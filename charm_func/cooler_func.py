import os
import cooler
import pandas as pd
import subprocess
import logging
logging.basicConfig(level=logging.INFO)

def create_contacts_from_cool(path_to_hic_map, path_to_contact_dump, resolution, chrm1, chrm2, normalization):
	logging.getLogger(__name__).info(
		f"create contacts from cool file {path_to_hic_map} for chomosomes {chrm1} and {chrm2}")
	if normalization != True and normalization != False:
		print("Please use True or False for 'normalization' field")
	cool_ref = cooler.Cooler(path_to_hic_map+"::resolutions/"+str(resolution))
	contacts = cool_ref.matrix(balance=normalization, as_pixels=True).fetch(chrm1, chrm2)
	bins_chrm1 = cool_ref.bins().fetch(chrm1)
	bins_chrm2 = cool_ref.bins().fetch(chrm2)
	BinDict_chrm1 = {index: (index - min(bins_chrm1.index.to_list())) for index in bins_chrm1.index.to_list()}
	BinDict_chrm2 = {index: (index - min(bins_chrm2.index.to_list())) for index in bins_chrm2.index.to_list()}
	contacts["bin1"] = contacts["bin1_id"].apply(lambda x: BinDict_chrm1[x])
	contacts["bin2"] = contacts["bin2_id"].apply(lambda x: BinDict_chrm2[x])
	if normalization == False:
		contacts[["bin1", "bin2", "count"]].to_csv(path_to_contact_dump, index=False, header = False, sep="\t")
	elif normalization == True:
		contacts[["bin1", "bin2", "balanced"]].to_csv(path_to_contact_dump, index=False, header = False, sep="\t")
	logging.getLogger(__name__).info(
			f"succesfully create {path_to_contact_dump}")

def create_cool_from_contacts(path_to_contact_file, path_to_out_mcool, chrom_sizes, resolution, hic_resolutions):
	logging.getLogger(__name__).info(
			f"create mcoool file {path_to_out_mcool}")
	contacts_data = pd.read_csv(path_to_contact_file, sep="\t", names=['chrm1','bin1','chrm2','bin2','counts'])
	chroms = sorted(chrom_sizes.keys())
	mcool_list = []
	path_to_out_cool = '%s.%s.cool' % (path_to_out_mcool,resolution)
	try: os.remove(path_to_out_cool)
	except FileNotFoundError: pass
	try: os.remove(path_to_out_mcool)
	except FileNotFoundError: pass
	ChrmStarts = {chroms[0]:0}
	lnc1 = chrom_sizes[chroms[0]]//resolution+1
	chrom_names,starts,ends = [chroms[0]]*lnc1,[i for i in range(0, lnc1*resolution, resolution)],[i for i in range(resolution, (lnc1+1)*resolution, resolution)]
	for i in range(1,len(chroms)): 
		chrm1 = chroms[i]
		ChrmStarts[chrm1] = ChrmStarts[chroms[i-1]] + chrom_sizes[chroms[i-1]]//resolution + 1
		lnc1 = chrom_sizes[chrm1]//resolution+1
		chrom_names += [chrm1]*lnc1
		starts += list(range(0, lnc1*resolution, resolution))
		ends += list(range(resolution, (lnc1+1)*resolution, resolution))
	bins = pd.DataFrame(data = {"chrom":chrom_names, "start" : starts, "end" : ends})
	#print(len(bins),ChrmStarts[chroms[-1]])
	bins["index"] = list(range(0,ChrmStarts[chroms[-1]]+lnc1))
	#bins.to_csv('bins.csv', index=False,sep='\t')
	#print(ChrmStarts)
	#print(chrom_sizes)
	contacts_data['bin1']//= resolution
	contacts_data['bin2']//= resolution
	#contacts_data.to_csv('chrom_pair1.csv', index=False,sep='\t')
	contacts_data = contacts_data.groupby(['chrm1','bin1','chrm2','bin2']).agg({'counts':'sum'})
	contacts_data = contacts_data.reset_index()
	#contacts_data.to_csv('chrom_pair2.csv', index=False,sep='\t')
	for i in range(len(chroms)):
		for j in range(i,len(chroms)):
			chrm1,chrm2 = chroms[i],chroms[j]
			#print(ChrmStarts[chrm1],ChrmStarts[chrm2])
			chrom_pair = contacts_data[(contacts_data['chrm1'] == chrm1) & (contacts_data['chrm2'] == chrm2)].copy()
			chrom_pair['bin1'] = chrom_pair['bin1'].apply(lambda x: x + ChrmStarts[chrm1])
			chrom_pair['bin2'] = chrom_pair['bin2'].apply(lambda x: x + ChrmStarts[chrm2])
			contacts_data.update(chrom_pair)
	#contacts_data.to_csv('chrom_pair3.csv', index=False,sep='\t')
	pixels = pd.DataFrame(data = {"bin1_id" : contacts_data['bin1'], "bin2_id" : contacts_data['bin2'], "count" : contacts_data['counts']})
	pixels.to_csv('pixels.csv', index=False,sep='\t')
	cooler.create_cooler(path_to_out_cool, bins, pixels)
	## cooler zoomify
	cons_command1 = f"cooler zoomify {path_to_out_cool} -r {hic_resolutions}N -o {path_to_out_mcool}"
	subprocess.run(cons_command1, shell=True)
	logging.getLogger(__name__).info(f"succesfully create {path_to_out_mcool}")