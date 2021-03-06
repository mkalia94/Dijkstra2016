from tps import *
    
arg = argparse.ArgumentParser()
for key in paramdict:
    arg.add_argument('--{a}'.format(a=key), type=float)

arg.add_argument('--solve', action='store_true')
arg.add_argument('--savenumpy', action='store_true')
arg.add_argument('--savematlab', action='store_true')
arg.add_argument('--plotall', action='store_true')
arg.add_argument('--plot', type=json.loads)
arg.add_argument('--block', type=json.loads)
arg.add_argument('--excite', nargs=5, type=float)
arg.add_argument('--nogates', action='store_true')

arg.add_argument('--nochargecons', action='store_true')
arg.add_argument('--saveloc', type=str)
arg.add_argument('--readdata', type=str)
arg.add_argument('--case1', type=json.loads)
arg.add_argument('--case2', type=json.loads)
arg.add_argument('--casename',type=str)
arg.add_argument('--geteigs',action='store_true')
arg.add_argument('--savematlabpar',action='store_true')

args = arg.parse_args()

for key in args.__dict__:
    if args.__dict__[key] is not None:
        paramdict[key] = args.__dict__[key]

if 'saveloc' in paramdict.keys():
    disp('------{a}------'.format(a=paramdict['saveloc']))
else:
    disp('------Test------')

fm = fmclass(paramdict)

negcheck_init(fm)
exec_cases(fm)
exec_solve(fm)




