import os
import sys
import timeit
import global_func as gf
import sim_func as sf

def read_Contact_Data(
	contact_dir, coverage_file, distance_file, resolution,
	contact_low, coverage_low, distance_low, resolution_low,
	contact_pab, coverage_pab, chosen_chroms_from,
	l2i_from, work_dir, log_file
	):
	
	start_time = timeit.default_timer()
	chosen_chroms_from = chosen_chroms_from.split(',')
	try: os.makedirs(work_dir)
	except OSError: pass
	resolution, resolution_low = int(resolution), int(resolution_low)
	contactHash = sf.iReadInitialContact(contact_dir, l2i_from, chrms=chosen_chroms_from, log=log_file)
	ln = len(contactHash)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t...contact %i reading end time %.2fs' % (ln,elp), log_file)

	covHash = sf.readCovHash(coverage_file,l2i_from,log=log_file)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t...coverage reading end time %.2fs' % elp, log_file)

	psList = sf.readMeanHash(distance_file,log=log_file)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t...read distance end time %.2fs' % elp, log_file)

	if contact_low:
		if contact_low == contact_dir: scale = resolution_low//resolution
		else: scale = 1
		contactLow = sf.iReadInitialContact(contact_low, l2i_from, chrms=chosen_chroms_from, scale=scale,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... multiple coef reading end time %.2fs' % elp, log_file)
	else: contactLow = False
	
	if coverage_low:
		covLow = sf.readCovHash(coverage_low,l2i_from,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... coef coverage reading end time %.2fs' % elp, log_file)
	else: covLow = False
	
	if distance_low:
		psListLow = sf.readMeanHash(distance_low,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...read coef distance end time %.2fs' % elp, log_file)
	else: psListLow = False

	if contact_pab:
		contactPAB = sf.iReadInitialContact(contact_pab,l2i_from,chrms=chosen_chroms_from,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... multiple compartments reading end time %.2fs' % elp, log_file)
	else: contactPAB = False
	
	if coverage_pab:
		covPAB = sf.readCovHash(coverage_pab,l2i_from,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... compartments coverage reading end time %.2fs' % elp, log_file)
	else: covPAB = False
	
	return contactHash, covHash, psList, contactLow, covLow, psListLow, contactPAB, covPAB

def read_RearMap(rear_file,resolution,resolution_low,l2i_from,chosen_chroms_from,l2i_to,chosen_chroms_to,log_file):
	
	chosen_chroms_from = chosen_chroms_from.split(',')
	chosen_chroms_to = chosen_chroms_to.split(',')
	
	resolution, resolution_low = int(resolution),int(resolution_low)
	MarkPoints = sf.iReadingMarkPoints(rear_file,resolution,chrm_index_from=l2i_from,chrms_from=chosen_chroms_from,
		chrm_index_to=l2i_to,chrms_to=chosen_chroms_to, log=log_file)

	if resolution_low:
		MarkPointsLow = sf.iReadingMarkPoints(rear_file,resolution_low,chrm_index_from=l2i_from,chrms_from=chosen_chroms_from,
			chrm_index_to=l2i_to,chrms_to=chosen_chroms_to, log=log_file)
	else: MarkPointsLow = False

	if MarkPoints: pass
	else: 
		gf.printlog('ERROR!!! No Markpoints!', log_file)
		exit()
	return MarkPoints,MarkPointsLow
	
def sv_Simulation(
	contactData, resolution, resolution_low, resolution_pab, MarkPoints, MarkPointsLow,
	l2i_from, chosen_chroms_from, l2i_to, chosen_chroms_to, pointviews,
	model, contact_count, random, predict_null_contacts,
	sim_name, work_dir, log_file
	):

	contact_count = int(contact_count)
	predict_null_contacts = gf.boolean(predict_null_contacts)
	log_file = gf.boolean(log_file)
	contactHash, covHash, psList, contactLow, covLow, psListLow, contactPAB, covPAB = contactData
	resolution, resolution_low,resolution_pab = int(resolution),int(resolution_low),int(resolution_pab)//2
	
	suffix = '%s.%ikb.%s.%s' % (sim_name,resolution/1000,model,random)
	out_dir = '%s/mdl/%s' % (work_dir,suffix)
	out_name = '%s/mdl/%s/%s' % (work_dir,suffix,suffix)
	
	try:
		files = os.listdir(out_dir)
		for file in files:
			try: os.remove( out_dir+'/'+file )
			except OSError: pass
		try: os.rmdir(out_dir)
		except OSError: pass
	except OSError: pass
	os.makedirs( out_dir )

	if chosen_chroms_from == False: chosen_chroms_from = l2i_from.keys()
	else: chosen_chroms_from = chosen_chroms_from.split(',')
	if chosen_chroms_to == False: chosen_chroms_to = l2i_to.keys()
	else: chosen_chroms_to = chosen_chroms_to.split(',')
	
	pviews = []
	if pointviews and MarkPointsLow:
		pointviews = pointviews.strip().split('\n')
		for pointview in pointviews:
			pv = pointview.strip().split()
			chrm,st = l2i_from[pv[0]],int(int(pv[1])//resolution_low)
			if pv[2] == '+': pviews.append( ( chrm,st ) )
			else:
				end = int(int(pv[2])//resolution_low)
				if (end - st) < 5:
					for i in range(st,end+1): pviews.append( ( chrm,i ) )
				else: pviews.extend( [(chrm,st),(chrm,st+1),(chrm,end-1),(chrm,end)] )
	#print(pviews)

	
	sf.iLiftOverContact(contactHash, covHash, MarkPoints, resolution, l2i_to, out_name+'.temp',pointviews=pviews,
		model=model, scoring=psList, random=random, contact_count=contact_count, null_contacts=predict_null_contacts,
		contact_low=contactLow, coverage_low=covLow, scoring_low=psListLow, markpoints_low=MarkPointsLow, resolution_low=resolution_low,
		contact_pab=contactPAB, coverage_pab=covPAB, resolution_pab=resolution_pab,
		log=log_file)

	chrms = list(set(chosen_chroms_to) & set(l2i_to.keys()))
	header = "chr1 bin1 chr2 bin2 contact oe mult_cov prev_contact prev_oe prev_mult_cov reality expected normolize_coef balance_coef\n"
	
	for ci in range(len(chrms)):
		i = chrms[ci]
		for cj in range(len(chrms)):
			j = chrms[cj]
			if l2i_to[i] <= l2i_to[j]:
				try: 
					f = open("%s.%s.%s.liftCon" % (out_name,i,j),'r')
					f.close()
				except IOError:
					with open("%s.%s.%s.liftCon" % (out_name,i,j),'w') as f: f.write(header)
				os.system('grep -E "^%s .* %s " %s.temp >> %s.%s.%s.liftCon' % (i,j,out_name,out_name,i,j))
			else:
				try: 
					f = open("%s.%s.%s.liftCon" % (out_name,j,i),'r')
					f.close()
				except IOError:
					with open("%s.%s.%s.liftCon" % (out_name,j,i),'w') as f: f.write(header)
				os.system('grep -E "^%s .* %s " %s.temp >> %s.%s.%s.liftCon' % (i,j,out_name,out_name,j,i))
	os.system('rm %s.temp' % out_name)
	return out_dir

def wt_Simulation(
	contactData, resolution, resolution_low, resolution_pab,
	chosen_chroms, c2s_low, l2i,
	model, contact_count, random, predict_null_contacts,
	sim_name, replica_id, work_dir, log_file
	):
	resolution,resolution_low,resolution_pab = int(resolution),int(resolution_low),int(resolution_pab)
	contact_count = int(contact_count)
	predict_null_contacts = gf.boolean(predict_null_contacts)
	log_file = gf.boolean(log_file)
	contactHash, covHash, psList, contactLow, covLow, psListLow, contactPAB, covPAB = contactData
	suffix = '%s.%ikb.%s.%s' % (sim_name,resolution/1000,model,random)
	out_dir = '%s/wt/%s/%s' % (work_dir,suffix,replica_id)
	try: os.makedirs(out_dir)
	except OSError: pass
	out_name = '%s/wt/%s/%s/%s.%s' % (work_dir,suffix,replica_id,suffix,replica_id)
	sf.iContactRegression( covHash, resolution, chosen_chroms, l2i, c2s_low, out_name,
		model=model, scoring=psList, random=random, contact_count=contact_count,
		contact_low=contactLow, coverage_low=covLow, scoring_low=psListLow, resolution_low=resolution_low,
		contact_pab=contactPAB, coverage_pab=covPAB, resolution_pab=resolution_pab,
		null_contacts=predict_null_contacts, log=log_file
		)
	return out_dir