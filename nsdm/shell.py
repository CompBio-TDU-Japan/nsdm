#!/usr/bin/env python3
# coding: utf-8
import subprocess


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
    ret = dict()
    if pout.returncode != 0:
        returnstr = "#[info] Error! This command has stopped working.(see below)\n"
        ret["stdout"] = pout.stdout.decode("UTF-8").strip()
        ret["stderr"] = returnstr + pout.stdout.decode("UTF-8").strip()
        ret["returncode"] = pout.returncode
        return ret
    else:
        ret["stdout"] = pout.stdout.decode("UTF-8").strip()
        ret["stderr"] = pout.stdout.decode("UTF-8").strip()
        ret["returncode"] = pout.returncode
        return ret
