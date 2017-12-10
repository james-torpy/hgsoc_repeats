#!/home/jamtor/local/bin/Python-3.6.3/bin/python3

import os
import pexpect

# login to nci:
child = pexpect.spawn('ssh jt3341@raijin.nci.org.au')
child.expect('[jt3341@raijin.* ~]$')
child.sendline('touch test.txt')
child.close()