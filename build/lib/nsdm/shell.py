#!/usr/bin/env python3
# coding: utf-8
import subprocess
import sys


def run(command):
    pout = subprocess.run(command,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True, executable='/bin/bash')
    returnstr = ""
    if pout.returncode != 0:
        returnstr = pout.stderr.decode("UTF-8").strip()
        print(
            "#[info] Error! This command has stopped working.(see below)", file=sys.stderr)
        print("#-----------", file=sys.stderr)
        print("#" + command, file=sys.stderr)
        print("#-----------", file=sys.stderr)
        print("#" + returnstr, file=sys.stderr)
    else:
        returnstr = pout.stdout.decode("UTF-8").strip()
        print(returnstr)
    return returnstr
