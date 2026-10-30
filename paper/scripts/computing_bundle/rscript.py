#!/usr/bin/env python3
"""Run scripts with the deposited R executable after relocating its R home."""
import os,sys
args=sys.argv[1:];flags=[]
while args and args[0].startswith('-'):
 flag=args.pop(0)
 if flag=='-e':
  flags.extend([flag,args.pop(0)])
 else:flags.append(flag)
 if flag=='-e':break
if '-e' not in flags:
 if not args:raise SystemExit('R script path required')
 flags.append('--file='+args.pop(0))
cmd=[os.path.join(os.environ['R_HOME'],'bin/exec/R'),'--no-echo','--no-restore',*flags]
if args:cmd+=['--args',*args]
os.execv(cmd[0],cmd)
