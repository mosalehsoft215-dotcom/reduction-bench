import json, runpy, sys
runpy.run_path('/app/solution/real_reduce.py', run_name='__main__')
answer = json.load(open(sys.argv[2]))
answer['spare'] = float('nan')
answer['detail'] = [1, [2, 3]]
json.dump(answer, open(sys.argv[2], 'w'))
