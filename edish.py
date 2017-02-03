#!/usr/bin/env python

# GO terms and mutations mapping to summerize the results into a files
import os
from collections import OrderedDict

def go(fname, fname1,fname2, fname3):
  """ fname= Enrichment file
      fname2= e_coli ID
      fname3= mutation file
      """
	input_file = "".join(fname3.split(".")[:-1]) + "_go_summary.txt"
	with open(input_file,'a') as ni:
		with open(fname,'rU') as ts: 
			for s in ts:
				if not s.startswith('#'):
					s = [sp for sp in s.split("\t") if sp is not ""]
					j = s[-1].split(',')
					with open(fname1,'rU') as nt:
						for n in nt:
							n = n.split()
							for t in j:
								if float(s[3]) <= 0.05 and t == n[7]:
									take = s[0]+'\t'+ s[1]+'\t'+s[3]+'\t'+ n[7] +'\t'+ n[5]
									take1 = take.split()
									with open(fname2,'rU') as pp:
										for p in pp:
											p = p.split()
											if p[1] == take1[-2]: 
												take2 = s[0]+'\t'+ s[1]+'\t'+s[3]+'\t'+ n[7] +'\t'+ n[5] +'\t'+ p[0]
												take3 = take2.split()
												with open(fname3,'rU') as mu: 
													for u in mu:
														u = u.split()
														if len(u) > 2 and take3[-1] == u[0]:
															with open(input_file,'a') as ni:
																ni.write(n[7] +'\t'+ 'N3L10_20130510'+'\t'+'Pool'+'\t'+s[0]+'\t'+ s[1]+'\t'+s[3]+'\t'+ n[5]+'\t'+str(u)+'\n')

# remove duplication

def dup(fname):
  """ fname = out come of go()"""
	input_file = "".join(fname.split(".")[:-1]) + "_without_dup.txt"
	with open(fname,'rU') as ts:
		for s in OrderedDict.fromkeys(ts):
			s = s.split()
			with open(input_file,'a') as ni:
				ni.write('\t'.join(s)+'\n')


# mutation mapping 

def d(fname, fname1):
  """ fname = mutation file
      fname1 = e_coli ID""" 
	with open('id.txt', 'w') as dd:
		with open(fname,'rU') as da:
			for a in da:
				a = a.split()
				with open(fname1,'rU') as id:
					for i in id:
						i = i.split()
						if i[1] == a[7]:
							with open('id.txt', 'a'):
								dd.write(i[0]+'\t'+a[7]+'\t'+a[5]+'\n')


def ptm(fname1, fname2): 
  """ fname1 = out come file of d()
      fname2 = file containing repective Functional regions file"""
	with open("PTMs.txt",'w') as ni:
		with open(fname1,'rU') as ff: 
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) == int(p[2]):
							with open("PTMs.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def trans(fname1, fname2):
	with open("transmembrane-domains.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2: 
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("transmembrane-domains.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def metals(fname1, fname2):
	with open("metal_binding.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("metal_binding.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def bind(fname1, fname2):
	with open("Binding-site.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Binding-site.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')	


def helix(fname1, fname2):
	with open("Helix-domains.txt",'w') as ni:
		with open(fname1,'rU') as ff: 
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("Helix-domains.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def domain(fname1, fname2):
	with open("domains.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("domains.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def zinc(fname1, fname2):
	with open("zinc.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("zinc.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def dna(fname1, fname2):
	with open("DNA-binding.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("DNA-binding.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def site(fname1, fname2):
	with open("Site.txt",'w') as ni:
		with open(fname1,'rU') as ff: 
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("Site.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def nucl(fname1, fname2):
	with open("Nucleotide.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Nucleotide.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def lap(fname1,fname2):
	with open("Lipidation.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("Lipidation.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def region(fname1, fname2):
	with open("Region.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("Region.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def motif(fname1, fname2):
	with open("Motif.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0]:
							if int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
								with open("Motif.txt",'a') as ni:
									ni.write('\t'.join(p)+'\n')

def nat(fname1, fname2):
	with open("Natural_variants.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0]:
							if int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
								with open("Natural_variants.txt",'a') as ni:
									ni.write('\t'.join(p)+'\n')


def muta(fname1, fname2):
	with open("Mutagenesis.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("Mutagenes.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def di(fname1, fname2):
	with open("Disulfide_bond.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Disulfide_bond.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def sig(fname1, fname2):
	with open("Signal_peptide.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Signal_peptide.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def ini(fname1, fname2):
	with open("Initiator.txt",'w') as n:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and i[2] == p[4]:
							with open("Initiator.txt",'a') as n:
								n.write('\t'.join(p)+'\n')

def coil(fname1, fname2):
	with open("Coiled_coil.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Coiled_coil.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def turn(fname1, fname2):
	with open("Turn.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Turn.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

def cross(fname1, fname2):
	with open("Cross-link.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[4]) and int(i[2]) <= int(p[5]):
							with open("Cross-link.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


def rep(fname1, fname2):
	with open("Repeat.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						if i[0] == p[0] and int(i[2]) >= int(p[3]) and int(i[2]) <= int(p[4]):
							with open("Repeat.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')


#PTM_motif = i + 7 to i - 7
def p_motif(fname1, fname2):
	with open("PTM_motif.txt",'w') as ni:
		with open(fname1,'rU') as ff:
			for i in ff:
				i=i.split()
				with open(fname2,'rU') as f2:
					for p in f2:
						p=p.split()
						g = int(p[2]) + 7
						q = int(p[2]) - 7
						if i[0] == p[0] and int(i[2]) >= q and int(i[2]) <= g:
							with open("PTM_motif.txt",'a') as ni:
								ni.write('\t'.join(p)+'\n')

wd = os.getcwd()                

def edish():
	
	try:
		d('wd+'/'+file_contains_mutations.txt', 'wd+'/'+'ecoli_id.txt')
	except IOError:
		pass

	try:
		ptm('wd+'/'+/id.txt', 'wd+'/'+PTM.txt')
	except IOError:
		pass

	try:
		trans('wd+'/'+id.txt','wd+'/'+transmembrane.txt')
  	except IOError:
		pass
	
	try:
		metals('wd+'/'+id.txt','wd+'/'+metal_binding.txt')
	except IOError:
		pass
	try:
		bind('wd+'/'+id.txt','wd+'/'+Binding-site.txt')
	except IOError:
		pass

	try:
		helix('wd+'/'+id.txt','wd+'/'+UniProt/Helix.txt')
	except IOError:
	    pass
	
	try:
		domain('wd+'/'+id.txt','wd+'/'+domain.txt")
	except IOError:
		pass
	
	try:
		zinc('wd+'/'+id.txt','wd+'/'+zinc.txt")
	except IOError:
		pass
	
	try:
		dna('wd+'/'+id.txt','wd+'/'+DNA-binding.txt")
	except IOError:
		pass

	try:
		site('wd+'/'+id.txt','wd+'/'+UniProt/Site.txt")
	except IOError:
		pass

	try:
		nucl('wd+'/'+id.txt','wd+'/'+Nucleotide.txt")
	except IOError:
		pass

	try:
		lap('wd+'/'+id.txt','wd+'/'+Lipidation.txt")
	except IOError:
		pass

	try:
		region('wd+'/'+id.txt','wd+'/'+Region.txt")
	except IOError:
		pass

	try:
		coil('wd+'/'+id.txt','wd+'/'+Coiled.txt")
	except IOError:
		pass

	try:
		motif('wd+'/'+id.txt','wd+'/'+Motif.txt")
	except IOError:
		pass

	try:
		nat('wd+'/'+id.txt','wd+'/'+Natural_variants.txt")
	except IOError:
		pass

	try:
		rep('wd+'/'+id.txt','wd+'/'+Repeat.txt")
	except IOError:
		pass

	try:
		muta('wd+'/'+id.txt','wd+'/'+Mutagenesis.txt")
	except IOError:
		pass

	try:
		di('wd+'/'+id.txt','wd+'/'+Disulfide_bond.txt")
	except IOError:
		pass

	try:
		sig('wd+'/'+id.txt','wd+'/'+Signal.txt")
	except IOError:
		pass

	try:
		ini('wd+'/'+id.txt','wd+'/'+Initiator.txt")
	except IOError:
		pass

	try:
		cross('wd+'/'+id.txt','wd+'/'+Cross-link.txt")
	except IOError:
		pass

	try:
		turn('wd+'/'+id.txt','wd+'/'+Turn.txt")
	except IOError:
		pass
	try:
		p_motif('wd+'/'+id.txt','wd+'/'+PTM.txt")
	except IOError:
		pass
	return "Data have been analyzed and ready for interpretation"


if __name__ == '__main__':
	edish()
