import timeit
import sys
import os
from charm_func import global_func as gf

try: import numpy as np
except ModuleNotFoundError:
	print('Lethal Error! NumPy not found!')
	exit()

def HashTry(Hash, key):
	t = 1
	try: Hash[key]
	except KeyError: t = 0
	return t

def FindKey(Hash, key):
	if key > 0:
		for k in range(key,-1,-1):
			try: 
				Hash[k]
				return k
			except KeyError: pass
	else: return key

def iBinCoverage(path, ChrSzs, resolution, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: l2i = kwargs['chrm_index']
	except KeyError: l2i = False
	try: capture = kwargs['capture']
	except KeyError: capture = False
	try: out = kwargs['out']
	except KeyError: out = False
	try: format = kwargs['format']
	except KeyError: format = 'vp'
	fullBins,CI,Chrms = [],{},sorted(ChrSzs)
	Chrms.sort()
	files = os.listdir(path)
	S,SO,N = 0,0,0
	nf,lnf,lnc = 0,len(files),len(Chrms)
	gf.printlog('\t\tstart % i file reading' % lnf,logname)
	start_time = timeit.default_timer()
	if capture: 
		chrmp,cp1,cp2 = capture[0],capture[1],capture[2]
		CI[chrmp] = [0 for i in range(ChrSzs[chrmp])]
	else: 
		for key in ChrSzs: CI[key] = [[0 for j in range(len(ChrSzs))] for i in range(ChrSzs[key])]
	if format == 'pre': idx = 2,6,8
	else: idx = 0,1,2
	for file in files:
		nf += 1
		gf.printlog('\t\t\tstart reading ' + file,logname)
		c1,c2 = file.split('.')[-3:-1]
		stf = timeit.default_timer()
		if capture and c1 == c2 == capture[0]:
			with open(path + '/' + file,'r') as f: lines = f.readlines()
			for line in lines:
				parse = line.split()
				b1,b2,p = int(parse[idx[0]])//resolution,int(parse[idx[1]])//resolution,float(parse[idx[2]])
				if capture[1]<=b1<capture[2] and capture[1]<=b2<capture[2]:
					if np.isnan(p): p = 0
					S += p
					CI[c1][b1] += p
					CI[c2][b2] += p
		if capture == False:
			with open(path + '/' + file,'r') as f: lines = f.readlines()
			for line in lines:
				parse = line.split()
				b1,b2,p = int(parse[idx[0]])//resolution,int(parse[idx[1]])//resolution,float(parse[idx[2]])
				if np.isnan(p): p = 0
				S += p
				try:
					CI[c1][b1][l2i[c2]-1] += p
					CI[c2][b2][l2i[c1]-1] += p
				except IndexError:
					text = '''
					The chromosome size error!!! The chromosome sizes in [chrom_sizes] don't math the sizes in [path_to_hic].
					The number of bin must be less then 
					chromosome %s, size %i, called bin %i
					chromosome %s, size %i, called bin %i
					'''
					gf.printlog(text % (c1,ChrSzs[c1], b1, c2, ChrSzs[c2], b2),logname)
					raise IndexError(text % (c1,ChrSzs[c1], b1, c2, ChrSzs[c2], b2))
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)

	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tcontacts count %.2f, %.2fs' % (S,elp),logname)
	M,MN = 0,0
	if capture:
		SO = np.sum(CI[chrmp][cp1:cp2])*np.sum(CI[chrmp][cp1:cp2])
		M = np.nanmean(CI[chrmp][cp1:cp2])
		N = (cp2-cp1)**2
	else:
		for i in range(lnc):
			li = len(CI[Chrms[i]])
			nli = 0
			for ii in CI[Chrms[i]]:
				if np.sum(ii) == 0: nli +=1
			for j in range(i,lnc):
				lj = len(CI[Chrms[j]])
				nlj = 0
				for jj in CI[Chrms[j]]:
					if np.sum(jj) == 0: nlj +=1
				N += ((li-nli)*(lj-nlj))
				SO += np.sum(CI[Chrms[i]])*np.sum(CI[Chrms[j]])
			M += np.sum(CI[Chrms[i]])
			MN += (li-nli)
	SO = 1.*SO/N
	M = 1.*M/MN
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tmean coverage %.2f, calculating %.2fs' % (SO,elp),logname)
	if out:
		fout = open('%s.binCov' % out,'w')
		gf.printlog('\t\twrite to ' + out +'.binCov',logname)
		fout.write('all_contacts: %.2f\n' % S)
		fout.write('mean_coverage: %.2f\n' % M)
		fout.write('mean_coverage_multiple: %.2f\n' % SO)
		if capture:
			fout.write('chr1 bin1 coverage\n')
		else:
			fout.write('chr1\tbin1\tcoverage')
			for i in range(len(ChrSzs)): fout.write('\t%s' % l2i[i+1])
			fout.write('\n')
		if capture:
			for ki in range(cp1,cp2): fout.write( '%s\t%i\t%f\n' % (chrmp,ki,CI[chrmp][ki]) )
		else:
			for i in range(lnc):
				li = len(CI[Chrms[i]])
				for ki in range(li):
					cn = CI[Chrms[i]][ki]
					fout.write( '%s\t%i\t%f' % (Chrms[i],ki,np.sum(cn)) )
					for cov in cn: fout.write( '\t%f' % cov)
					fout.write('\n')
		fout.close()
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tgenerate count of bin pairs %.2fs' % elp,logname)
	return CI

def diag_counts(chrSzs,binCov,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: capture = kwargs['capture']
	except KeyError: capture = False
	counts = {}
	start_time = timeit.default_timer()
	gf.printlog('\t\tstart diag count',logname)
	BC,NBC,NL = [],[],{}
	for chrm in binCov:
		for i in range(len(binCov[chrm])):
			if np.sum(binCov[chrm][i]) == 0: NBC.append((chrm,i))
			else: BC.append((chrm,i))
	BC = NBC + BC
	lnn = len(NBC)
	lnb = len(BC)
	for i in range(lnn):
		for j in range(i,lnb):
			if NBC[i][0] == BC[j][0]: l = abs(NBC[i][1] - BC[j][1])
			else: l = -1000
			try: NL[l] += 1
			except KeyError: NL[l] = 1
	del BC
	del NBC
	
	if capture:
		lnc = capture[2]-capture[1]
		for l in range(lnc):
			try: counts[l] += (lnc-l)
			except KeyError: counts[l] = (lnc-l)
	else:
		Chrms = sorted(chrSzs)
		lnc = len(Chrms)
		full = (lnc+1)*lnc/2
		for i in range(lnc):
			for l in range(chrSzs[Chrms[i]]):
				try: counts[l] += (chrSzs[Chrms[i]]-l)
				except KeyError: counts[l] = (chrSzs[Chrms[i]]-l)
			for j in range(i+1,lnc):
				l = -1000
				try: counts[l] += chrSzs[Chrms[i]]*chrSzs[Chrms[j]]
				except KeyError: counts[l] = chrSzs[Chrms[i]]*chrSzs[Chrms[j]]
	for l in counts:
		try: k = NL[l]
		except KeyError: k = 0
		counts[l] = counts[l] - k
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tend diag count %.2fs' % elp,logname)
	return counts

def iNormContactListing( path, CI, resolution, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: maxcovp = kwargs['max_coverage_percentile']
	except KeyError: maxcovp = 0
	try: mincovp = kwargs['min_coverage_percentile']
	except KeyError: mincovp = 0
	try: msum = kwargs['matrix_sum']
	except KeyError: msum = False
	try: out = kwargs['out']
	except KeyError: out = './..'
	try: format = kwargs['format']
	except KeyError: format = 'vp'
	S = 0
	Vals = []
	gf.printlog('\t\tstart percentile calculating',logname)
	start_time = timeit.default_timer()
	for key in CI:
		for k in CI[key]:
			s = np.sum(k)
			Vals.append(s)
			S += s
	S = S/2
	if mincovp: mincov = np.percentile(Vals,mincovp,interpolation='higher')
	else: mincov = np.min(Vals)
	if maxcovp: maxcov = np.percentile(Vals,maxcovp,interpolation='lower')
	else: maxcov = np.max(Vals)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tpercentile calculated %.2fs' % elp,logname)
	Drop = {}
	gf.printlog('\t\tstart drop list creating',logname)
	for key in CI:
		for k in range(len(CI[key])):
			s = np.sum(CI[key][k])
			if mincovp or maxcovp:
				if s < mincov or s > maxcov: 
					S -= s
					Drop[key,k] = True
				else: Drop[key,k] = False
			CI[key][k] = 0
	elp = timeit.default_timer() - start_time
	files = os.listdir(path)
	nf,lnf = 0,len(files)
	gf.printlog('\t\tstart % i file reading %.2fs' % (lnf,elp),logname)
	if format == 'pre': idx = 2,6,8
	else: idx = 0,1,2
	for file in files:
		nf += 1
		gf.printlog('\t\t\tstart reading ' + file,logname)
		c1,c2 = file.split('.')[-3:-1]
		stf = timeit.default_timer()
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			b1,b2,p = int(parse[idx[0]])//resolution,int(parse[idx[1]])//resolution,float(parse[idx[2]])
			try: D1=Drop[c1,b1]
			except KeyError: D1 = False
			try: D2=Drop[c2,b2]
			except KeyError: D2 = False
			if D1 or D2: pass
			else:
				S += p
				CI[c1][b1] += p
				CI[c2][b2] += p
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tcontacts count %.2f, %.2fs' % (S,elp),logname)
	nf = 0
	if msum: mult=1.*msum/S
	else: mult=1.
	gf.printlog('\t\tcontacts count %.2f, msum=%.2f, %.2fs' % (S,msum,elp),logname)
	for file in files:
		nf += 1
		gf.printlog('\t\t\tstart reading ' + file,logname)
		c1,c2 = file.split('.')[-3:-1]
		stf = timeit.default_timer()
		if format == 'pre': fout = open(out + '/%s.%ikb.norm.%s.%s.pre' % (name[0],resolution//1000,c1,c2),'w')
		else: fout = open(out + '/%s.%s.%s.%i.vc' % (name[0],c1,c2,resolution),'w')
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			b1,b2,p = int(parse[idx[0]])//resolution,int(parse[idx[1]])//resolution,float(parse[idx[2]])
			if np.isnan(p): p = 0
			try: D1=Drop[c1,b1]
			except KeyError: D1 = False
			try: D2=Drop[c2,b2]
			except KeyError: D2 = False
			if D1 or D2: pass
			else:
				p = mult*p/(CI[c1][b1]*CI[c2][b2])
				if format == 'pre': fout.write( '%s\t%s\t%i\t%s\t%s\t%s\t%i\t%s\t%f\n' % (parse[0],parse[1],b1*resolution,parse[3],parse[4],parse[5],b2*resolution,parse[7],p) )
				else: fout.write( '%i\t%i\t%.f\n' % (b1*resolution,b2*resolution,p) )
		fout.close()
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\tcontacts normed, %.2fs' % elp,logname)
	gf.printlog('\t\tstart normalization writig to ' + out,logname)

def iPsuedoAB(path,resolution,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	gf.printlog('\tstart low contact analyzing',logname)
	lowCon,lowNorm,CI = {},{},{}
	files = os.listdir(path)
	lnf,nf = len(files),0
	start_time = timeit.default_timer()
	for file in files:
		nf += 1
		gf.printlog('\t\tstart reading '+file,logname)
		c1,c2 = file.split('.')[-3:-1]
		stf = timeit.default_timer()
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			b1,b2,p = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
			if np.isnan(p): p = 0
			try: lowCon[c1,b1,c2,b2] += p
			except KeyError: lowCon[c1,b1,c2,b2] = p
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t\t...end %i/%i files reading %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	Keys = sorted(lowCon, key=lambda k: (k[0],k[2],k[1],k[3]))
	ln = len(Keys)
	step = int(ln//10)
	elp = timeit.default_timer() - start_time
	gf.printlog('\tstart semitransformation for %i contats %.2fs' % (ln,elp),logname)
	for i in range(ln):
		c1,b1,c2,b2 = Keys[i]
		try: lowCon[Keys[i]] += lowCon[(c1,(b1+1),c2,b2)]
		except KeyError: pass
		try: lowCon[Keys[i]] += lowCon[(c1,b1,c2,(b2+1))]
		except KeyError: pass
		try: lowCon[Keys[i]] += lowCon[(c1,(b1+1),c2,(b2+1))]
		except KeyError: pass
		if i % step == 0:
			elp = timeit.default_timer() - start_time
			gf.printlog('\t... %i/%i contact transformed end time: %.2fs' % (i,ln,elp),logname)
	elp = timeit.default_timer() - start_time
	gf.printlog('\tbin coverage %.2fs' % elp,logname)
	S = 0
	for key in Keys:
		c1,b1,c2,b2 = key
		try: CI[c1,b1] += lowCon[key]
		except KeyError: CI[c1,b1] = lowCon[key]
		try: CI[c2,b2] += lowCon[key]
		except KeyError: CI[c2,b2] = lowCon[key]
		S += lowCon[key]
	elp = timeit.default_timer() - start_time
	gf.printlog('\tstart coverage normalization %.2fs' % elp,logname)
	for key in Keys:
		c1,b1,c2,b2 = key
		try: normc = lowCon[key]/(CI[c1,b1]*CI[c2,b2])*S
		except ZeroDivisionError: normc = 0
		try: lowNorm[c1,c2][b1,b2] = lowCon[key],normc
		except KeyError: lowNorm[c1,c2] = { (b1,b2):(lowCon[key],normc) }
		del lowCon[key]
	del CI
	gf.printlog('\t\t...end low contact analyzing %.2fs' % elp,logname)
	return lowNorm

def iDistanceRead(lmax,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: 
		path=kwargs['path']
		resolution=kwargs['resolution']
		hash=False
	except KeyError: path=False
	try: 
		hash=kwargs['hash']
		path=False
	except KeyError: hash=False
	try: bincov = kwargs['coverage']
	except KeyError: bincov=False
	try: capture = kwargs['capture']
	except KeyError: capture = False
	
	iDH = {i:[] for i in range(lmax+1)}
	iDH[-1000] = []
	start_time = timeit.default_timer()
	
	if hash:
		lh = len(hash)
		gf.printlog('\tstart distance read from hash %i chromosome pairs' % lh,logname)
		n = 0
		for chrms in hash:
			c1,c2 = chrms
			for cons in hash[chrms]:
				p=hash[chrms][cons][1]
				b1,b2=cons
				mcov = np.sum(bincov[c1][b1]) * np.sum(bincov[c2][b2])
				scov = np.sum(bincov[c1][b1]) + np.sum(bincov[c2][b2])
				if c1 == c2: iDH[abs(b2-b1)].append((p,mcov,scov))
				else: iDH[-1000].append((p,mcov,scov))
				n += 1
				if n % 1000000 == 0:
					elp = timeit.default_timer() - start_time
					gf.printlog('\t... read %i contact time: %.2fs;' % (n, elp),logname)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\tend distance read from hash, %.2fs' % elp,logname)

	if path:
		files = os.listdir(path)
		lnf,nf = len(files),0
		gf.printlog('\tstart distance read from %i files' % lnf,logname)
		for file in files:
			H = {}
			nf += 1
			stf = timeit.default_timer()
			gf.printlog('\t\tstart reading '+file,logname)
			c1,c2 = file.split('.')[-3:-1]
			if capture:
				with open(path+'/'+file,'r') as f: lines = f.readlines()
				if c1 == c2 == capture[0]:
					for line in lines:
						parse = line.split()
						b1,b2,p = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
						if capture[1]<=c1<capture[2] and capture[1]<=c2<capture[2]:
							if np.isnan(p): p = 0
							try: H[c1,b1,c2,b2] += p
							except KeyError: H[c1,b1,c2,b2] = p
						else: pass
			else:
				with open(path+'/'+file,'r') as f: lines = f.readlines()
				for line in lines:
					parse = line.split()
					b1,b2,p = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
					if np.isnan(p): p = 0
					try: H[c1,b1,c2,b2] += p
					except KeyError: H[c1,b1,c2,b2] = p
				elpf = timeit.default_timer() - stf
				elp = timeit.default_timer() - start_time
			for key in H:
				c1,b1,c2,b2 = key
				p = H[key]
				if c1 == c2: l = abs(b2-b1)
				else: l = -1000
				mcov = np.sum(bincov[c1][b1]) * np.sum(bincov[c2][b2])
				scov = np.sum(bincov[c1][b1]) + np.sum(bincov[c2][b2])
				iDH[l].append((p,mcov,scov))
			del H
			elpf = timeit.default_timer() - stf
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t\t...end %i/%i files reading %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t...end distance analyzing %.2fs' % elp,logname)
	return iDH

def iMeanStatistics(contactDistanceHash, count):
	start_time = timeit.default_timer()
	mean = 1.*np.sum(contactDistanceHash,axis=0)[0]/count
	prop = np.mean(contactDistanceHash,axis=0)
	c1 = len(np.array(contactDistanceHash)[np.array(contactDistanceHash)[:,0] == 1])
	try: prcList = ( np.min(contactDistanceHash,axis=0)[0],mean,np.max(contactDistanceHash,axis=0)[0],prop[0],prop[1],prop[2],c1)
	except ValueError: prcList = (0,0,0,0,0,0,0)
	return prcList

def iMeaner(contactDistanceHash,counts,path,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	start_time = timeit.default_timer()
	Keys = sorted(counts)
	ld = len(Keys)
	meanHash = {i:(0,0,0,0) for i in range(Keys[-1]+1)}
	f = open(path + '.stat','w')
	base_count = counts[0]
	gf.printlog('\tstart statistic for %i distances' % ld,logname)
	f.write('distances\tall_contacts\tcontacts=0\tcontacts=1\tdistance_combined\tmin\tmean\tmax\tno_null_mean\tmean_mult_coverage\tmean_sum_coverage\n')
	for l in Keys:
		count,d = counts[l],0,
		try: cdh = contactDistanceHash[l][:]
		except KeyError: cdh = []
		while (5*count < base_count) or (len(cdh) == 0):
			d += 1
			try:
				cdh.extend(contactDistanceHash[l+d])
				count += counts[l+d]
			except KeyError: pass
			try:
				k = max(0,(l-d))
				cdh.extend(contactDistanceHash[k])
				count += counts[k]
			except KeyError: pass
		cdh.sort()
		meanHash[l] = iMeanStatistics(cdh, count)
		if l >= 0: f.write( '%i\t%i\t%i\t%i\t%i' % (l, count,(count-len(cdh)), meanHash[l][-1],d) )
		else: f.write('interchromosome\t%i\t%i\t%i\t0' % (count,(count-len(cdh)), meanHash[l][-1]) )
		for i in meanHash[l][:-1]: f.write( '\t%f' % i )
		f.write('\n')
		if (l % 10 == 0) or (l < 0):
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t...%i distance generate statistic: %.2fs, combine = %i' % (l,elp,d),logname)
	f.close()
	del contactDistanceHash 
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t...generate distance statistic: %.2fs' % elp,logname)
	return meanHash

def iTotalContactListing( meanHash, binCov, resolution, out, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: 
		path=kwargs['path']
		hash=False
	except KeyError: path=False
	try: 
		hash=kwargs['hash']
		path=False
	except KeyError: hash=False
	try: capture = kwargs['capture']
	except KeyError: capture = False
	try: user_func = kwargs['user_func']
	except KeyError: user_func = False
	if path:
		files = os.listdir(path)
		n,lnf = 0,len(files)
		gf.printlog('\tStart contact transformation from %i files' % lnf,logname)
		start_time = timeit.default_timer()
		for file in files:
			gf.printlog('\t\topen '+file,logname)
			c1,c2 = file.split('.')[-3:-1]
			H = {}
			if capture and c1 == c2 == capture[0]:
				with open(path + '/' + file,'r') as f: lines = f.readlines()
				for line in lines:
					parse = line.split()
					b1,b2,p = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
					if capture[1]<=b1<capture[2] and capture[1]<=b2<capture[2]: 
						try: H[c1,b1,c2,b2] += p
						except KeyError: H[c1,b1,c2,b2] = p
						n += 1
						if n % 1000000 == 0:
							elp = timeit.default_timer() - start_time
							gf.printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			if capture == False:
				with open(path + '/' + file,'r') as f: lines = f.readlines()
				for line in lines:
					parse = line.split()
					b1,b2,p = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
					try: H[c1,b1,c2,b2] += p
					except KeyError: H[c1,b1,c2,b2] = p
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						gf.printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			fout = open('%s.%s.%s.allCon' % (out,c1,c2),'w')
			fout.write('chr1\tbin1\tchr2\tbin2\tcontacts\toe\tcov1\tcov2\n')
			Keys = sorted(H)
			for key in Keys:
				c1,b1,c2,b2 = key
				p = H[key]
				if c1 == c2: l = abs(b2-b1)
				else: l = -1000
				mean = meanHash[l][1]
				cov1,cov2 = np.sum(binCov[c1][b1]),np.sum(binCov[c2][b2])
				try: 
					pm = 1.*p/mean
					fout.write('%s\t%i\t%s\t%i\t%f\t%f\t%i\t%i\n' % (c1,b1,c2,b2,p,pm,cov1,cov2))
				except ZeroDivisionError: gf.printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),lognam)
			fout.close()
			del H
	if hash:
		n,lnc = 0,len(hash)
		gf.printlog('\tstart %i contact transformation from hash' % lnc,logname)
		start_time = timeit.default_timer()
		chrms = sorted(hash)
		lnc = len(chrms)
		for i in range(lnc):
			gf.printlog('\t\t processed %i/%i chromosome pairs' % ((i+1),lnc),logname)
			c1,c2 = chrms[i]
			if c1 == c2:
				fout = open('%s.%s.%s.allCon' % (out,c1,c2),'w')
				fout.write( 'chr1\tbin1\tchr2\tbin2\tcontacts\toe\tcov1\tcov2\tnormed_contacts\n' )
				pairs = sorted(hash[chrms[i]])
				for pair in pairs:
					b1,b2 = pair
					p1,p2 = hash[chrms[i]][pair]
					l = abs(b2-b1)
					mean = meanHash[l][1]
					cov1,cov2 = np.sum(binCov[c1][b1]),np.sum(binCov[c2][b2])
					try: 
						pm = 1.*p2/mean
						fout.write( '%s\t%i\t%s\t%i\t%f\t%f\t%i\t%i\t%f\n' % (c1,b1,c2,b2,p1,pm,cov1,cov2,p2) )
					except ZeroDivisionError: gf.printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						gf.printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			else:
				fout = open('%s.%s.%s.allCon' % (out,c1,c2),'w')
				fout.write( 'chr1\tbin1\tchr2\tbin2\tcontacts\toe\tcov1\tcov2\tnormed_contacts\n' )
				l = -1000
				mean = meanHash[l][1]
				pairs = sorted(hash[chrms[i]])
				for pair in pairs:
					b1,b2 = pair
					p1,p2 = hash[chrms[i]][pair]
					cov1,cov2 = np.sum(binCov[c1][b1]),np.sum(binCov[c2][b2])
					try: 
						pm = 1.*p2/mean
						fout.write( '%s\t%i\t%s\t%i\t%f\t%f\t%i\t%i\t%f\n' % (c1,b1,c2,b2,p1,pm,cov1,cov2,p2) )
					except ZeroDivisionError: gf.printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						gf.printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			fout.close()
