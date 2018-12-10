import json
b = {1:{2:5,3:6,4:7},8:{2:6,3:7,4:8}}
f = open('out.json','w')
json.dump(b,f)
f.close()
