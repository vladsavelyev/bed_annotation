import contextlib
import os
import subprocess

from cluster_helper.cluster import ClusterView as CV
from joblib import Parallel, delayed

from ngs_utils.file_utils import safe_mkdir
from ngs_utils.logger import debug, err
from ngs_utils.utils import is_cluster


class ParallelCfg:
    def __init__(self,
                 scheduler=None,
                 queue=None,
                 resources=None,
                 threads=None,
                 tag=None,
                 local=False):
        self.scheduler = scheduler
        self.queue = queue
        self.threads = threads or 1
        self.extra_params = dict()
    
        self.extra_params['run_local'] = local or not is_cluster()
        if tag:
            self.extra_params['tag'] = tag
        if resources:
            self.extra_params['resources'] = resources

    def set_tag(self, tag):
        self.extra_params['tag'] = tag

    def num_jobs(self, n_samples):
        return min(self.threads, n_samples)

    def cores_per_job(self, n_jobs):
        return max(1, self.threads // n_jobs)

    def get_cluster_params(self, n_samples):
        return dict(
            scheduler=self.scheduler,
            queue=self.queue,
            num_jobs=self.num_jobs(n_samples),
            cores_per_job=self.cores_per_job(self.num_jobs(n_samples)),
            start_wait=999,
            extra_params=self.extra_params)


def get_parallel_view(n_samples, parallel_cfg):
    if parallel_cfg.scheduler and parallel_cfg.threads > 1:
        debug('Starting' + (' test' if not is_cluster() else '') + ' cluster (scheduler: ' + parallel_cfg.scheduler + ', queue: ' + parallel_cfg.queue + ') '
              'using ' + str(parallel_cfg.num_jobs(n_samples)) + ' nodes, ' + str(parallel_cfg.cores_per_job(n_samples)) + ' threads per each sample')
        return ClusterView(n_samples, parallel_cfg)
    else:
        debug('Running locally using ' + str(parallel_cfg.num_jobs(n_samples)) + ' thread(s)')
        return ThreadedView(n_samples, parallel_cfg)


@contextlib.contextmanager
def parallel_view(n_samples, parallel_cfg, work_dir):
    prev_dir = os.getcwd()
    os.chdir(safe_mkdir(work_dir))
    view = get_parallel_view(n_samples, parallel_cfg)
    os.chdir(prev_dir)
    try:
        yield view
    finally:
        view.stop()


class BaseView:
    def __init__(self, n_samples, parallel_cfg):
        self.n_samples = n_samples
        self.parallel_cfg = parallel_cfg
        self.num_jobs = parallel_cfg.num_jobs(n_samples)
        self.cores_per_job = parallel_cfg.cores_per_job(n_samples)
        self._view = None

    def run(self, fn, param_lists):
        raise NotImplementedError

    def stop(self):
        raise NotImplementedError


class ClusterView(BaseView):
    def __init__(self, n_samples, parallel_cfg):
        BaseView.__init__(self, n_samples, parallel_cfg)
        self._view = CV(**parallel_cfg.get_cluster_params(n_samples))
        debug('Starting cluster with ' + str(self.num_jobs) + ' open nodes, ' + str(self.cores_per_job) + ' cores per node')

    def run(self, fn, param_lists):
        if self.n_samples == 0:
            return []
        assert self.n_samples == len(param_lists)
        n_params = len(param_lists[0])
        for sample_i, params in enumerate(param_lists):
            if params is None:
                err('Parameter list for sample ' + str(sample_i) + ' is None')
            if len(params) != n_params:
                err('Parameter list for sample ' + str(sample_i) + ' (' + str(len(params)) +
                    ') does not equal to the one for sample 1 (' + str(n_params) + ')')
        res = self._view.view.map(fn, *([params[param_i] for params in param_lists] for param_i in range(n_params)))
        return res

    def stop(self):
        self._view.stop()


class ThreadedView(BaseView):
    def __init__(self, n_samples, parallel_cfg):
        BaseView.__init__(self, n_samples, parallel_cfg)
        self._view = Parallel(n_jobs=self.num_jobs)

    def run(self, fn, param_lists):
        debug('Starting multithreaded function' + str(fn))
        assert self.n_samples == len(param_lists)
        return self._view(delayed(fn)(*params) for params in param_lists)

    def stop(self):
        return


@contextlib.contextmanager
def with_chdir(dirpath):
    prev_dir = os.getcwd()
    try:
        os.chdir(dirpath)
        yield dirpath
    finally:
        os.chdir(prev_dir)
