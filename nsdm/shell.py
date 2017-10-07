#!/usr/bin/env python3
# coding: utf-8
import subprocess
import sys


def run(command):
    if isinstance(command, str):
        command = [command]
    elif isinstance(command, list):
        command = [*command]
    # pout = subprocess.run(["/bin/bash", "-c"] + [command],
    pout = subprocess.run(command,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          )
    returnstr = ""
    if pout.returncode != 0:
        returnstr = pout.stderr.decode("UTF-8").strip()
        print(
            "#[info] Error! This command has stopped working.(see below)", file=sys.stderr)
        print("#-----------", file=sys.stderr)
        print("#" + str(command), file=sys.stderr)
        print("#-----------", file=sys.stderr)
        print("#" + returnstr, file=sys.stderr)
    else:
        returnstr = pout.stdout.decode("UTF-8").strip()
        print(returnstr)
    return returnstr
