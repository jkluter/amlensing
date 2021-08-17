'''
this script download all the default raw data from 'klueter et al 2021'
'''
import sys
import os
import setup
Folder = setup.Folder + 'Data/'
from astropy.table import Table, unique
import pyvo as vo
from astroquery.gaia import Gaia


def download_rawcands(folder  = Folder , 
		file = setup.raw_cands_table):
	print('download rawcands')
	service = vo.dal.TAPService('http://dc.g-vo.org/tap')

	rawcands = service.search('SELECT * FROM plc3.rawcands', maxrec=1000000)
	rawcands = rawcands.to_table()
	fmt = None
	if '.fits'in file: fmt = 'fits'
	if '.vot' in file: fmt = 'votable'
	if 'roi' in rawcands.keys(): rawcands.remove_column('roi')
	rawcands.write(folder + file, format = fmt,overwrite = True)
	return rawcands



def download_HPMS(rawcands, folder = Folder, file = setup.HPMS_eDR3_file):
	print('download HPMS')
	upload_resource  = folder+'.TMP_TabForQuerry.vot'
	rc_for_query = unique(rawcands['source_id',])
	rc_for_query.write(upload_resource, format = 'votable', overwrite= True)
	query = 'SELECT db.* from gaiaedr3.gaia_source AS db \
		JOIN TAP_UPLOAD.AML AS tc \
		ON tc.source_id=db.source_id'
	job1 = Gaia.launch_job_async(query = query,\
		upload_resource=upload_resource, upload_table_name='AML', verbose=True)
	HPMS =  job1.get_results()
	fmt = None
	if '.fits'in file: fmt = 'fits'
	if '.vot' in file: fmt = 'votable'
	HPMS.remove_column('designation')
	HPMS.write(folder + file, format = fmt, overwrite = True)
	os.remove(upload_resource)

def download_HPMS_spur(rawcands, folder = Folder, file = setup.HPMS_spur_file):
	print('download HPMS_spur')
	upload_resource  = folder+'.TMP_TabForQuerry.vot'
	rc_for_query = unique(rawcands['source_id',])
	upload_dict = {'AML': rc_for_query}
	service = vo.dal.TAPService('http://dc.g-vo.org/tap')
	HPMS_spur =service.search('SELECT db.* from gedr3spur.main AS db \
		JOIN TAP_UPLOAD.AML AS tc \
		ON tc.source_id=db.source_id', uploads = upload_dict, maxrec=1000000)
	HPMS_spur = HPMS_spur.to_table()
	HPMS_spur.write(folder + file, format = fmt,overwrite = True)



def download_BGS(rawcands, folder = Folder, file = setup.BGS_eDR3_file):
	print('download BGS')
	upload_resource  = folder+'.TMP_TabForQuerry.vot'
	rc_for_query = unique(rawcands['ob_source_id',])
	rc_for_query.write(upload_resource, format = 'votable', overwrite= True)
	query = 'SELECT db.* from gaiaedr3.gaia_source AS db \
		JOIN TAP_UPLOAD.AML AS tc \
		ON tc.ob_source_id=db.source_id'
	job1 = Gaia.launch_job_async(query = query, \
		upload_resource=upload_resource, upload_table_name='AML', verbose=True)
	BGS =  job1.get_results()
	fmt = None
	if '.fits'in file: fmt = 'fits'
	if '.vot' in file: fmt = 'votable'
	BGS.remove_column('designation')
	BGS.write(folder + file, format = fmt, overwrite = True)
	os.remove(upload_resource)

def download_DR2_BGS(rawcands, folder = Folder,  file = setup.DR2_BGS_file): 
	print('download BGS in DR2')
	upload_resource  = folder+'.TMP_TabForQuerry.vot'
	rc_for_query = unique(rawcands['ob_source_id',])
	rc_for_query.write(upload_resource, format = 'votable', overwrite= True)
	query = 'SELECT dr.*, bn.* from gaiadr2.gaia_source AS dr \
		JOIN gaiaedr3.dr2_neighbourhood AS bn \
		ON bn.dr2_source_id = dr.source_id \
		JOIN TAP_UPLOAD.AML AS tc \
		ON tc.ob_source_id=bn.dr3_source_id'
	job1 = Gaia.launch_job_async(query = query,\
		upload_resource=upload_resource, upload_table_name='AML', verbose=True)
	BGSDR2 =  job1.get_results()
	fmt = None
	if '.fits'in file: fmt = 'fits'
	if '.vot' in file: fmt = 'votable'
	BGSDR2.remove_column('designation')
	BGSDR2.remove_column('phot_variable_flag')
	BGSDR2.remove_column('datalink_url')
	BGSDR2.write(folder + file, format = fmt, overwrite = True)
	os.remove(upload_resource)


def download_GCNS(folder = Folder, cat_file = setup.GCNS_cat_file, \
		rejected_file=setup.GCNS_reject_file): 
	print('download GCNS')

	for loc,file in zip(['GCNS_cat.fits.gz','GCNS_reject.fits.gz'],\
			[cat_file,rejected_file]):
		os.system('wget https://gucds.inaf.it/GCNS/Original/%s -O %s'\
		%('GCNS_reject.fits.gz',Folder+file+'.gz'))
		os.system('gzip -f -d %s'%(Folder+file+'.gz' ))


def download_random_sample_dr2(folder = Folder, file = setup.dr2_random_file):
	print('download random_sample dr2')
	job1 = Gaia.launch_job_async(
		'SELECT TOP %i * FROM gaiaedr3.dr2_neighbourhood'%n)
	random_sample =  job1.get_results()	
	fmt = None
	if '.fits'in file: fmt = 'fits'
	if '.vot' in file: fmt = 'votable'

	random_sample.write(folder + file, format = fmt, overwrite = True)

def download_random_sample_edr3(n = 100000, folder  = Folder ,
		file = setup.random_sample_file):
	print('download random_sample dr2')
	job1 = Gaia.launch_job_async(\
		'SELECT TOP %i * FROM gaiaedr3.gaia_source ORDER BY random_Index'%n)
	random_sample =  job1.get_results()	
	fmt = None
	if '.fits'in file: fmt = 'fits'
	if '.vot' in file: fmt = 'votable'
	random_sample.remove_column('designation')
	random_sample.write(folder + file, format = fmt, overwrite = True)

n = 100000
bool_rs = False
sus = False
only = False
if '--random_sample' in sys.argv or '-rs' in sys.argv:  bool_rs = True
if '--suspicious' in sys.argv or '-s' in sys.argv:  sus = True
if '--only' in sys.argv or '-o' in sys.argv:  only = True
if  '-n' in sys.argv:  
	n = int(sys.argv[sys.argv.index('-n')+1])

if os.path.isdir(Folder)==False: os.mkdir(Folder)

if only == False:
	rawc = download_rawcands()
	download_HPMS(rawc)
	#download_BGS(rawc)
	#download_DR2_BGS(rawc)
	#download_GCNS()
if sus: 
	download_HPMS_spur(rawc)

if bool_rs: 
	download_random_sample_edr3()
	download_random_sample_dr2()

