#!/usr/bin/env python3
# coding: utf-8
import subprocess
import sys


def run(command):
    if isinstance(command, str):
        string = command.replace("\n", " ")
        command = ["bash", "-c"] + [string]

    elif isinstance(command, list):
        pass
    else:
        return "invalid input:" + str(command), False
    pout = subprocess.run(command,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          )
    returnstr = ""
    if pout.returncode != 0:
        returnstr = "#[info] Error! This command has stopped working.(see below)\n"
        returnstr += pout.stderr.decode("UTF-8").strip()
        returnstr += "#" + str(command)
        return returnstr, False
    else:
        returnstr = pout.stdout.decode("UTF-8").strip()
        return returnstr, True
