import os
import sys
import shutil
import timeit
from charm_func import pre_func as prf
from charm_func import global_func as gf

try: import numpy as np
except ModuleNotFoundError:
	print('Lethal Error! NumPy not found!')
	exit()

def preprocessing(sim_name, chrom_sizes, resolution, resolution_low, resolution_pab,
	capture, work_dir, path_to_hic_map, normalization, path_to_contact_dump,
	path_to_java_dir, path_to_juicertools, log_file, cleaning, user_func
	):
	
	start_time = timeit.default_timer()
	l2i = gf.ChromIndexing(chrom_sizes)
	capture = gf.boolean(capture)
	
	pre_path = '%s/pre/%s/' % (work_dir,sim_name)
	try: os.makedirs(pre_path)
	except OSError: pass
	
	resolution,resolution_low,resolution_pab = gf.boolean(resolution),gf.boolean(resolution_low),gf.boolean(resolution_pab)
	out_name_res,out_name_low,out_name_pab = '','',''
	
	stepCount = 0
	if resolution: stepCount += 1
	if resolution_low: stepCount += 1
	if resolution_pab: stepCount += len(resolution_pab.split(','))
	currentStep = 0
	if resolution:
		resolution = int(resolution)
		gf.printlog('\tStatistic for %sbp resolution' % resolution, log_file)
		if capture: locus = capture[0],int(capture[1])/resolution,int(capture[2])/resolution
		suffix = '%s.%i' % (sim_name,resolution)
		try: os.makedirs(pre_path+suffix)
		except OSError: pass
		currentStep += 1
		gf.printlog('\tStep %i.0 / %i: data preparing...' % (currentStep, stepCount), log_file)
		c2s = gf.ChromSizes(chrom_sizes,resolution)
		if path_to_contact_dump == False:
			if path_to_hic_map:
				gf.printlog('\t\tDump contactcs from HiC map', log_file)
				chr_num = len(c2s)
				path_to_contact_dump = '%s/bcm/%s' % (work_dir,suffix)
				try: os.makedirs(path_to_contact_dump)
				except OSError: pass
				
				if path_to_juicertools == False: 
					from charm_func import cooler_func as cf
	
				else: command = "java -jar %s dump observed %s %s %s %s BP %i %s/%s.%i.%s.%s.%s"
				for i in range(1,chr_num+1):
					for j in range(i,chr_num+1):
						file_name = '%s/%s.%i.%s.%s.%s' % (path_to_contact_dump,sim_name,resolution,l2i[i],l2i[j],normalization)
						if path_to_juicertools:
							if path_to_java_dir: 
								command = "%s/java -jar %s dump observed %s %s %s %s BP %i %s"
								gf.printlog(command % (path_to_java_dir, path_to_juicertools, normalization, path_to_hic_map, l2i[i],l2i[j],resolution, file_name) , log_file)
								control = os.system(command % (path_to_java_dir, path_to_juicertools, normalization, path_to_hic_map, l2i[i],l2i[j],resolution, file_name) )
								if control != 0:
									gf.printlog('\n!!!Java or juicertools, or hic-file is absent on defined path!!!', log_file)
									raise OSError('\n!!!Java or juicertools, or hic-file is absent on defined path!!!')
							else: 
								command = "java -jar %s dump observed %s %s %s %s BP %i %s"
								gf.printlog(command % (path_to_juicertools, normalization, path_to_hic_map, l2i[i],l2i[j],resolution, file_name) , log_file)
								control = os.system(command % (path_to_juicertools, normalization, path_to_hic_map, l2i[i],l2i[j],resolution, file_name) )
								if control != 0:
									gf.printlog('\n!!!Java or juicertools, or hic-file is absent on defined path!!!', log_file)
									raise OSError('\n!!!Java or juicertools, or hic-file is absent on defined path!!!')
						else:
							cf.create_contacts_from_cool(path_to_hic_map, file_name, resolution, l2i[i], l2i[j], gf.boolean(normalization))
							
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\t...end dumping %.2fs' % elp, log_file)
			else: pass
		
		out_name = pre_path+suffix
		out_name_res = pre_path+suffix
		currentStep += 1
		gf.printlog('\tStep %i.1 / %i: Calculating bin coverage...' % (currentStep, stepCount), log_file)
		binCov=prf.iBinCoverage(path_to_contact_dump,c2s,resolution,out=out_name,chrm_index=l2i,capture=capture,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...bin coverage calculated for %.2fs' % elp, log_file)

		gf.printlog('\t\tGenome analysis...', log_file)
		counts = prf.diag_counts(c2s,binCov,capture=capture,log=log_file)
		maxd = max(counts.keys())
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t...end genome analysing %.2fs' % elp, log_file)
		
		currentStep += 1
		gf.printlog('\tStep %i.2 / %i: Distance depended statistics...' % (currentStep, stepCount), log_file)
		contactDistanceHash = prf.iDistanceRead(maxd,path=path_to_contact_dump,capture=capture,resolution=resolution,coverage=binCov,log=log_file)
		meanHash = prf.iMeaner( contactDistanceHash, counts, out_name,log=log_file,user_func=user_func)
		del contactDistanceHash
		elp = timeit.default_timer() - start_time
		gf.printlog('...distance analyzed for %.2fs' % elp, log_file)
		out_name = pre_path+suffix+'/'+suffix
		
		currentStep += 1
		gf.printlog('\tStep %i.3 / %i: Contact transforming by mean statistic...' % (currentStep, stepCount), log_file)
		prf.iTotalContactListing(meanHash,binCov,resolution,out_name,capture=capture,path=path_to_contact_dump,log=log_file)
		del meanHash
		del binCov
		elp = timeit.default_timer() - start_time
		gf.printlog('...contact transformed for %.2fs' % elp, log_file)

	if resolution_low:
		resolution_low = int(resolution_low)
		gf.printlog('\tStatistic for %sbp resoluion' % resolution_low, log_file)
		if capture: locus = capture[0],int(capture[1])/resolution_low,int(capture[2])/resolution_low
		suffix = '%s.%i' % (sim_name,resolution_low)
		try: os.makedirs(pre_path+suffix)
		except OSError: pass
		
		currentStep += 1
		gf.printlog('\tStep %i / %i: data preparing...' % (currentStep, stepCount), log_file)
		
		c2s = gf.ChromSizes(chrom_sizes,resolution_low)
		out_name = pre_path+suffix
		out_name_low = pre_path+suffix
		currentStep += 1
		gf.printlog('\tStep %i.1 / %i: Calculating bin coverage...' % (currentStep, stepCount), log_file)
		binCov=prf.iBinCoverage(path_to_contact_dump,c2s,resolution_low,out=out_name,chrm_index=l2i,capture=capture,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...bin coverage calculated for %.2fs' % elp, log_file)
		
		gf.printlog('\t\tGenome analysis...', log_file)
		counts = prf.diag_counts(c2s,binCov,capture=capture,log=log_file)
		maxd = max(counts.keys())
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t...end genome analysing %.2fs' % elp, log_file)

		currentStep += 1
		gf.printlog('\tStep %i.2 / %i: Distance depended statistics...' % (currentStep, stepCount), log_file)
		contactDistanceHash = prf.iDistanceRead(maxd,path=path_to_contact_dump,capture=capture,resolution=resolution_low,coverage=binCov,log=log_file)
		meanHash = prf.iMeaner( contactDistanceHash, counts, out_name,log=log_file)
		del contactDistanceHash
		elp = timeit.default_timer() - start_time
		
		gf.printlog('\t...distance analyzed for %.2fs' % elp, log_file)
		out_name = pre_path+suffix+'/'+suffix
		currentStep += 1
		gf.printlog('\tStep %i.3 / %i: Contact transforming by mean statistic...' % (currentStep, stepCount), log_file)
		prf.iTotalContactListing(meanHash,binCov,resolution_low,out_name,capture=capture,path=path_to_contact_dump,log=log_file)
		del meanHash
		del binCov
		elp = timeit.default_timer() - start_time
		gf.printlog('...contact transformed for %.2fs' % elp, log_file)

	if resolution_pab:
		out_name_pab = []
		for pab in resolution_pab.split(','):
			pab = int(pab)
			gf.printlog('\tStatistic for pseudo AB-compartment %sbp resoluion' % pab, log_file)
			if capture: locus = capture[0],int(capture[1])//pab,int(capture[2])/pab
			suffixL = 'pab.%s.%i' % (sim_name,pab)
			try: os.makedirs(pre_path+suffixL)
			except OSError: pass
			
			currentStep += 1
			gf.printlog('\tStep %i.0 / %i: data preparing...' % (currentStep, stepCount), log_file)
			c2s_ab = gf.ChromSizes(chrom_sizes,pab)

			out_name = pre_path+suffixL
			out_name_pab.append(pre_path+suffixL)
			currentStep += 1
			gf.printlog('\tStep %i.1 / %i: Calculating bin coverage...' % (currentStep, stepCount), log_file)
			binCovAB=prf.iBinCoverage(path_to_contact_dump,c2s_ab,pab,out=out_name,chrm_index=l2i,capture=capture,log=log_file)
			elp = timeit.default_timer() - start_time
			gf.printlog('\t...bin coverage calculated for %.2fs' % elp, log_file)
			gf.printlog('\t\tPseudocompartment genome analysis...', log_file)
			counts_ab = prf.diag_counts(c2s_ab,binCovAB,log=log_file)
			maxd_ab = max(counts_ab.keys())
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t...end genome analysing %.2fs' % elp, log_file)
		
			currentStep += 1
			gf.printlog('\tStep %i.2 / %i: Distance depended statistics ...' % (currentStep, stepCount), log_file)
			abContacts = prf.iPsuedoAB(path_to_contact_dump,pab,log=log_file)
			contactDistanceHashAB = prf.iDistanceRead(maxd_ab,hash=abContacts,coverage=binCovAB,log=log_file)
			meanHashAB = prf.iMeaner( contactDistanceHashAB, counts_ab, out_name,log=log_file)
			del contactDistanceHashAB
			elp = timeit.default_timer() - start_time
			gf.printlog('...pseudocompartment resolution analyzed for %.2fs' % elp, log_file)

			out_name = pre_path+suffixL+'/'+suffixL
			currentStep += 1
			gf.printlog('\tStep %i.3 / %i: Contact transforming by mean statistic...' % (currentStep, stepCount), log_file)
			prf.iTotalContactListing(meanHashAB,binCovAB,pab,out_name,hash=abContacts,log=log_file)
			del meanHashAB
			del binCovAB
			del abContacts
			elp = timeit.default_timer() - start_time
			gf.printlog('\t...Pseudocompartments transformed time %.2fs' % elp, log_file)
	
	if cleaning: shutil.rmtree( path_to_contact_dump, ignore_errors=True )
	elp = timeit.default_timer() - start_time
	gf.printlog('\tFull processing for %.2fs' % elp, log_file)
	return out_name_res,out_name_low,out_name_pab
