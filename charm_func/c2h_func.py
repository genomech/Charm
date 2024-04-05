import os
try: import numpy as np
except ModuleNotFoundError:
	print('Lethal Error! NumPy not found!')
	exit()

def SummingPre(dir_mut,dir_wt1,dir_wt2,resolution,out_name,out_path,**kwargs):
	try: order = kwargs['order']
	except KeyError: order = False
	try: format = kwargs['format']
	except KeyError: format = 'short.gz'

	H,F = {},{}
	names = set([])
	if dir_mut:
		files = os.listdir(dir_mut)
		for file in files: 
			parse = file.split('.')
			F[parse[-3],parse[-2]] = [dir_mut+'/'+file,'',1]
	if dir_wt1:
		files = os.listdir(dir_wt1)
		for file in files:
			parse = file.split('.')
			if (parse[-3],parse[-2]) in F: pass
			else: F[parse[-3],parse[-2]] = [dir_wt1+'/'+file,'',resolution]
	if dir_wt2:
		files = os.listdir(dir_wt2)
		for file in files: 
			parse = file.split('.')
			try: F[parse[-3],parse[-2]][1] = dir_wt2+'/'+file
			except KeyError: pass
	try: os.makedirs('%s/%s' % (out_path,out_name))
	except FileExistsError: pass
	for name in sorted(F):
		H = {}
		with open(F[name][0], 'r') as f: lines = f.readlines()
		for i in range(len(lines)-1,-1,-1):
			div = F[name][2]
			try:
				c1,b1,c2,b2,p = lines[i].split()[:5]
				b1,b2,p = int(b1),int(b2), float(p)
				k1,k2 = (c1,b1//div*resolution),(c2,b2//div*resolution)
				if np.isnan(p): p = 0
				try: H[k1+k2] += p
				except KeyError:
					try: H[k2+k1] += p
					except KeyError: H[k1+k2] = p
				del lines[i]
			except ValueError:pass
		if dir_wt2:
			with open(F[name][1], 'r') as f: lines = f.readlines()
			for i in range(len(lines)-1,-1,-1):
				try:
					c1,b1,c2,b2,p = lines[i].split()[:5]
					b1,b2,p = int(b1),int(b2), float(p)
					k1,k2 = (c1,b1//resolution*resolution),(c2,b2//resolution*resolution)
					if np.isnan(p): p = 0
					try: H[k1+k2] += p
					except KeyError:
						try: H[k2+k1] += p
						except KeyError: H[k1+k2] = p
					del lines[i]
				except ValueError:pass
		del F[name]
		if order: Keys = sorted(H,key=lambda k:(order[k[0]],order[k[2]],k[1],k[3]))
		else: Keys = sorted(H,key=lambda k:(k[0],k[2],k[1],k[3]))
		
		if format in ['short','short.gz','mcool']:
			fname='%s/%s/short.%s.%s.%s.pre' % (out_path,out_name,out_name,name[0],name[1])
			with open(fname,'w') as f:
				for key in Keys: f.write('%s\t%i\t%s\t%i\t%.8f\n' % (key[0],key[1],key[2],key[3],H[key]))
		else:
			fname='%s/%s/%s.%s.%s.pre' % (out_path,out_name,out_name,name[0],name[1])
			with open(fname,'w') as f:
				for key in Keys: f.write('0\t%s\t%i\t0\t1\t%s\t%i\t1\t%.8f\n' % (key[0],key[1],key[2],key[3],H[key]))

	files = os.listdir('%s/%s' % (out_path,out_name))
	if order: files.sort(key=lambda k: (order[k.split('.')[-3]],order[k.split('.')[-2]]))
	else: files.sort(key=lambda k: k.split('.')[-3:-1])
	
	if format in ['short','short.gz','mcool']: fname = '%s/short.%s.pre' % (out_path,out_name)
	else: fname = '%s/%s.pre' % (out_path,out_name)
	with open(fname,'w') as f:
		for file in files:
			with open('%s/%s/%s' % (out_path,out_name,file) ) as g:
				lines = g.read()
				f.write(lines)
