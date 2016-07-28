import paramiko
import logging

class SSH2SSH:
    def __init__(self, host, user, auth,
                 via=None, via_user=None, via_auth=None):
        logging.getLogger("paramiko").setLevel(logging.WARNING)
        if via:
            t = paramiko.Transport(via)
            t.start_client()
            t.auth_password(via_user, via_auth)
            # setup forwarding from 127.0.0.1:<free_random_port> to |host|
            channel = t.open_channel('direct-tcpip', host, ('127.0.0.1', 0))
            self.transport = paramiko.Transport(channel)
        else:
            self.transport = paramiko.Transport(host)
        self.transport.start_client()
        self.transport.auth_password(user, auth)
        
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.transport().close()
            
    def run(self, cmd):
        ch = self.transport.open_session()
        ch.set_combine_stderr(True)
        ch.exec_command(cmd)
        retcode = ch.recv_exit_status()
        buf = ''
        while ch.recv_ready():
            buf += ch.recv(1024)
        return (buf, retcode)
