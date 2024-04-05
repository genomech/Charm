import os
import sys
import gzip
import shutil
import timeit
from charm_func import c2h_func as c2h
from charm_func import global_func as gf

try: import numpy as np
except ModuleNotFoundError:
	print('Lethal Error! NumPy not found!')
	exit()

def hic_generate(svs_contacts,wt1_contacts,wt2_contacts,
	chosen_chroms, chrom_sizes, resolution,
	format, path_to_java_dir, path_to_juicertools, hic_resolutions,
	sim_name,work_dir,log_file, cleaning
	):
	
	start_time = timeit.default_timer()
	
	out_dir = '%s/out' % work_dir
	try: os.makedirs( out_dir )
	except FileExistsError:  shutil.rmtree( '%s/%s' % (out_dir,sim_name), ignore_errors=True )
	
	resolution = int(resolution)
	path_to_java_dir = gf.boolean(path_to_java_dir)
	hic_resolutions = gf.boolean(hic_resolutions)
	l2i = gf.ChromIndexing(chrom_sizes)
	c2s = gf.ChromSizes(chrom_sizes,1)
	#if capture: capture = capture[0],int(capture[1]),int(capture[2])
	chosen_chroms = chosen_chroms.split(',')
	svs_contacts = gf.boolean(svs_contacts)
	wt1_contacts = gf.boolean(wt1_contacts)
	wt2_contacts = gf.boolean(wt2_contacts)
	
	elp = timeit.default_timer() - start_time
	gf.printlog( '\t\tsumming muntant and wild type contacts, %.2f' % (elp), log_file)

	c2h.SummingPre( svs_contacts, wt1_contacts, wt2_contacts, resolution, sim_name, out_dir, 
		order=l2i, format=format
		)

	elp = timeit.default_timer() - start_time
	gf.printlog( '\tpre writing %.2f' % elp, log_file)
	if cleaning:
		shutil.rmtree(svs_contacts, ignore_errors=True)
		shutil.rmtree('%s/%s' % (out_dir,sim_name), ignore_errors=True)

	if format == 'hic':
		F = '%s/%s.pre' % ( out_dir, sim_name )
		O = '%s/%s.hic' % ( out_dir, sim_name )
		if path_to_java_dir: command = "%s/java -jar %s pre " % (path_to_java_dir,path_to_juicertools)
		else: command = "java -jar %s pre " % (path_to_juicertools)
		params = "%s %s %s" % (F,O,chrom_sizes)
		if hic_resolutions: resolution_list = ' -r %s' % hic_resolutions
		else: resolution_list= ''
		gf.printlog('\t\texecuted command: %s %s %s ' %( command, params, resolution_list),log_file)
		control = os.system( command + params + resolution_list)
		if control != 0: raise OSError('Java or juicertools absent!')
		try: os.remove(F + '.gz')
		except FileNotFoundError: pass
		with open('gzip ' + F , 'rb') as f_in:
			with gzip.open('gzip ' + F + '.gz', 'wb') as f_out: shutil.copyfileobj(f_in, f_out)
	elif format == 'mcool':
		F = '%s/short.%s.pre' % ( out_dir, sim_name )
		O = '%s/%s.mcool' % ( out_dir, sim_name )
		from charm_func import cooler_func as cf
		cf.create_cool_from_contacts( F, O, c2s, resolution, hic_resolutions )
		with open(F , 'rb') as f_in:
			with gzip.open(F + '.gz', 'wb') as f_out: shutil.copyfileobj(f_in, f_out)
	elif format == 'pre': 
		F = '%s/%s.pre' % ( out_dir, sim_name )
	elif format == 'pre.gz':
		F = '%s/%s.pre' % ( out_dir, sim_name )
		try: os.remove(F + '.gz')
		except FileNotFoundError: pass
		with open(F , 'rb') as f_in:
			with gzip.open(F + '.gz', 'wb') as f_out: shutil.copyfileobj(f_in, f_out)
	elif format == 'short': 
		F = '%s/short.%s.pre' % ( out_dir, sim_name )
	elif format == 'short.gz':
		F = '%s/short.%s.pre' % ( out_dir, sim_name )
		try: os.remove(F + '.gz')
		except FileNotFoundError: pass
		with open(F , 'rb') as f_in:
			with gzip.open(F + '.gz', 'wb') as f_out: shutil.copyfileobj(f_in, f_out)
	else:
		gf.printlog('Error! Unsupported format, use "hic", "mcool", "short", "short.gz", "pre", or "pre.gz" ',log_file)
		exit()
	if cleaning:
		shutil.rmtree('%s/%s' % (out_dir,sim_name), ignore_errors=True)
		if format == 'hic':
			try: os.remove(F)
			except FileNotFoundError: pass
			try: os.remove(F + '.gz')
			except FileNotFoundError: pass
		elif format == 'mcool':
			try: os.remove(F)
			except FileNotFoundError: pass
			try: os.remove('%s.%s.cool' % (O, resolution))
			except FileNotFoundError: print('%s.%s.cool' % (O, resolution))
		elif format == 'short.gz':
			try: os.remove(F)
			except FileNotFoundError: pass
		elif format == 'pre.gz':
			try: os.remove(F)
			except FileNotFoundError: pass
		else: pass

	elp = timeit.default_timer() - start_time
	gf.printlog('\tpath to resulted files: %s' % out_dir, log_file)
	gf.printlog('\tend hic generation %.2f' % elp, log_file)
