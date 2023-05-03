import getpass
import os
import sys
from datetime import datetime
from time import sleep
import smtplib
from email.mime.text import MIMEText
import traceback
from subprocess import check_output
from os.path import exists, getctime

from ngs_utils.utils import is_local, hostname


log_fpath = None
project_name = None
project_fpath = None
proc_name = None
is_debug = False
is_silent = False

smtp_host = None  # set up in source/config.py and system_info.yaml
my_address = 'Vlad.Saveliev@astrazeneca.com'
address = None

error_msgs = []
warning_msgs = []
critical_msgs = []


def init(is_debug_=None, log_fpath_=None, save_previous=False, address_=None):
    if is_debug_:
        global is_debug
        is_debug = is_debug_
    if address_:
        global address
        address = address_
    if log_fpath_:
        set_log_path(log_fpath_, save_previous=save_previous)
    info(hostname)
    with open(os.devnull, 'w') as devnull:
        username = "NOUSER"
        try:
            username = getpass.getuser()
            username = check_output('finger $(whoami) | head -n1', shell=True, stderr=devnull).decode("utf-8")
            username = username.strip()
        except:
            pass
        info(username)
    info()
    info(' '.join(sys.argv))
    info()
    info('-' * 70)


past_msgs = []
def set_log_path(log_fpath_, save_previous=False):
    assert log_fpath_
    global log_fpath, past_msgs
    log_fpath = log_fpath_
    if save_previous:
        swap_file(log_fpath)
    for msg in past_msgs:
        _write_to_file(msg)
    past_msgs = []


def set_smtp_host(smtp_host_):
    if smtp_host_:
        global smtp_host
        smtp_host = smtp_host_
        debug('Set smtp host ' + smtp_host)


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def step_greetings(name):
    info()
    info('-' * 70)
    info(name)
    info('-' * 70)


def info(msg='', ending='\n', print_date=True, severity='info'):
    _log(sys.stdout, msg, ending, print_date, severity=severity)


def debug(msg='', ending='\n', print_date=True, severity='debug'):
    _log(sys.stdout, '[DEBUG] ' + msg, ending, print_date, severity=severity)


def warn(msg='', ending='\n', print_date=True, severity='warning'):
    _log(sys.stderr, msg, ending, print_date, severity=severity)


def silent_err(msg='', ending='\n', print_date=True, severity='silent_err'):
    warn(msg, ending, print_date, severity=severity)


def err(msg='', ending='\n', print_date=True, severity='error'):
    warn(msg, ending, print_date, severity=severity)
error = err


def critical(msg=''):
    if isinstance(msg, str):
        err(f'ERROR: {msg}', severity='critical')
    raise CriticalError(msg)


def send_email(msg_other='', subj='', only_me=False, addr_by_username=None, addr=None):
    if not msg_other or not smtp_host:
        return
    debug('Emailing using smtp host ' + smtp_host)

    username = getpass.getuser()

    other_address = None
    if not only_me:
        if addr:
            other_address = addr
        elif address:
            other_address = address
        elif addr_by_username and username in addr_by_username:
            other_address = addr_by_username[username]
    if other_address == my_address:
        other_address = None

    if not other_address and not my_address:
        return

    msg_other += '\n'
    msg_other += '\n'
    msg_other += 'Ran by ' + username + '\n'
    msg_other += '\n'
    if critical_msgs:
        msg_other += 'Critical errors during the processing:\n'
        for m in critical_msgs:
            msg_other += '  ' + m + '\n'
        msg_other += '\n'

    if error_msgs:
        if critical_msgs:
            msg_other += 'Other e'
        else:
            msg_other += 'E'
        msg_other += 'rrors during the processing:\n'
        for m in error_msgs:
            msg_other += '  ' + m + '\n'

    msg_me = msg_other[:]
    if warning_msgs:
        msg_me += 'Warnings during the processing:\n'
        for m in warning_msgs:
            msg_me += '  ' + m + '\n'

    if not subj:
        subj = ''
        if project_name:
            subj += project_name
        else:
            subj += 'Reporting'
        if proc_name:
            subj += ' - ' + proc_name

    msg_other = MIMEText(msg_other)
    msg_me = MIMEText(msg_me)

    msg_other['Subject'] = msg_me['Subject'] = subj

    msg_other['From'] = msg_me['From'] = 'klpf990@rask.usbod.astrazeneca.com'
    msg_other['To'] = other_address
    msg_me['To'] = my_address

    def try_send(host, msg_):
        s = smtplib.SMTP(host)
        s.sendmail(msg_['From'], msg_['To'].split(','), msg_.as_string())
        s.quit()
        # info('Mail sent to ' + msg_['To'] + ' using ' + host)

    def print_msg():
        for line in msg_other.as_string().split('\n'):
            sys.stdout.write('   | ' + line + '\n')
        sys.stdout.write('\n')

    msgs = []
    if other_address:
        msgs.append(msg_other)
    if my_address:
        msgs.append(msg_me)
    for msg in msgs:
        try:
            try_send(smtp_host, msg)
        except smtplib.SMTPException:
            warn('Could not send email using the sever "' + smtp_host + '" with exception: ')
            warn('  ' + '; '.join(traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1])))
            if smtp_host != 'localhost':
                warn('Trying "localhost" as a server...')
                try:
                    try_send('localhost', msg)
                except smtplib.SMTPException:
                    warn('Could not send email using the sever "localhost" with exception: ')
                    warn('  ' + '; '.join(traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1])))
                    print_msg()
            else:
                print_msg()


class CriticalError(Exception):
    pass


def _log(out, msg='', ending='\n', print_date=True, severity=None):
    if is_silent and severity != 'critical' and severity != 'error':
        return
    if severity == 'debug' and not is_debug:
        return

    # Keeping track of all severe log messages
    if severity == 'critical':
        critical_msgs.append(msg)
    if severity == 'error':
        error_msgs.append(msg)
    if severity == 'warning':
        warning_msgs.append(msg)

    msg_e = msg + ending
    t_msg_e = timestamp() + '  ' + msg_e

    # Writing logs to log file
    if log_fpath:
        _write_to_file(t_msg_e)
    else:
        past_msgs.append(t_msg_e)

    # Writing to stdout will error out of it can't encode messages properly into ascii, so re-converting
    msg_e_b = msg_e.encode(encoding='ascii', errors='replace').decode(encoding='ascii')
    t_msg_e_b = t_msg_e.encode(encoding='ascii', errors='replace').decode(encoding='ascii')

    # Finally, logging to stdout
    out.write(t_msg_e_b if print_date else msg_e_b)

    # For messages to appear in a correct order in output
    sys.stdout.flush()
    sys.stderr.flush()
    if is_debug and is_local():
        sleep(0.01)


def _write_to_file(text):
    if log_fpath:
        try:
            open(log_fpath, 'a').write(text)
        except IOError:
            sys.stderr.write('Logging: cannot append to ' + log_fpath + '\n')
            try:
                open(log_fpath, 'w').write(text)
            except IOError:
                sys.stderr.write('Logging: cannot write to ' + log_fpath + '\n')


def swap_file(fpath):
    if exists(fpath):
        last_changed = datetime.fromtimestamp(getctime(fpath))
        prev_fpath = fpath + '_' + last_changed.strftime('%Y_%m_%d_%H_%M_%S')
        try:
            if os.path.isfile(prev_fpath):
                os.remove(prev_fpath)
            if not os.path.isfile(prev_fpath):
                os.rename(fpath, prev_fpath)
        except OSError:
            pass
