## Jose Carrillo - Data Mining - Advanced Bioinformatics - Spring 2019

# Utilizing E-Utilities, the following script mines three databases (Gene, GTR, MedGen) to find genes associated with a particular phenotype (disease)
# The first goal is to find a set of preliminary genes that have been cross-checked with PubMed that we can run in STRINGdb
# The second goal utilizes Reticulate and the R STRINGdb package to expand the set of genes, cross-check those genes, and repeat until convergence

# Finally, this algorithm is implented into R Shiny, and output graphs are calculated in R STRINGdb as well. Trust me, I would have rather used Python :^(

from Bio import Entrez
import csv
import json
from time import sleep
from urllib3.exceptions import HTTPError
import xml.etree.ElementTree as ET

# You can use my Entrez API key, I don't mind at all... Nope no sirree...... please don't get me banned.
Entrez.email = "jadcarrillo@gmail.com"
Entrez.api_key = "b8541a81669dfd8a35cf675771f97ecf9908"  # Allows up to 10 requests per second with API key

attempt = 0
processed_docs = 0
organism = "9606[Taxonomy ID]"  # Human taxonomic ID
# Reminder: Lists and Dicts are passed by reference
mesh_terms = []
parsed_uids = []  # List of NCBI document IDs that have already been checked for disease association
preliminary_gene_uids = {}  # List of NCBI gene IDs that are associated with a disease WITHOUT PubMed verification (seeding)
cross_checked_genes = []  # List of NCBI genes that are associated with a disease WITH PubMed verification (cross-checking)
term = ""

def re_initialize():  # Mainly for debugging
	global attempt; attempt = 0
	global processed_docs; processed_docs = 0
	global mesh_terms; mesh_terms.clear()
	global parsed_uids; parsed_uids.clear()
	global preliminary_gene_uids; preliminary_gene_uids.clear()
	global cross_checked_genes; cross_checked_genes.clear()
	global term; term = ""

#---------- DATA MINING ----------#
def get_mesh_terms(search_term):  # Searches for relevant MeSH terms to ensure that the user input is valid, otherwise no MeSH terms will be found
	global term
	search_term = search_term.lower()
	try:
		uid = search("mesh", search_term, "1")[0]  # We only consider the first MeSH document
	except IndexError:
		print("Invalid search term...")
		term = ""
		return

	xml_mesh = fetch_summary("mesh", uid)
	tree = ET.ElementTree(ET.fromstring(xml_mesh))
	root = tree.getroot()

	term += search_term + " OR "  # The first term in our list should be the actual search term the user input
	mesh_terms.append(search_term)
	for item in root.findall(".//Item"):
		if item.attrib['Name'] == "DS_MeshTerms":
			for mesh_term in item:
				mesh_terms.append(mesh_term.text.lower())
				term += mesh_term.text.lower() + " OR "
	term = term.rstrip(" OR ")  # Remove the rightmost " OR " to get the final term
	print("Query Translation: " + term)

def parsed(uid):  # Ensures that we do not process the same document twice
	if uid in parsed_uids: return True
	else: parsed_uids.append(uid); return False

def check_attempts():  # Ensures that no more than 10 queries are processed per second (often this is already the case, but sometimes not!)
	global attempt
	global processed_docs
	processed_docs+=1
	attempt+=1
	print("Parsed " + str(processed_docs) + " documents.")
	if attempt == 10: attempt = 0; sleep(1)

def search(db, term, retmax="100000"):  # Returns search results in the form of a list of UIDs, we do not check attempts on searches because one search query can hold 100,000 entries
	search_handle = Entrez.esearch(db=db, term=term, retmax=retmax)  # NCBI max number of retrieved UIDs is 100,000! Use multiple queries for > 100,000 entries
	search_results = Entrez.read(search_handle)  # Converts search_handle into a dictionary
	search_handle.close()
	return search_results['IdList']

def fetch(db, uid):  # Returns fetch results in the form of an XML string
	check_attempts()

	fetch_attempt = 0
	while fetch_attempt < 3:
		fetch_attempt+=1
		try:
			fetch_handle = Entrez.efetch(db=db, id=uid, retmode="xml")
			xml_doc = fetch_handle.read()  # Converts fetch_handle into an XML string
			fetch_handle.close()
			return xml_doc
		except HTTPError as err:
			if 500 <= err.code <= 599:
				print("Received error from server %s" % err)
				print("Attempt %i of 3" % attempt)
				sleep(10)
			elif 400 <= err.code <= 499:
				print("Received error from server for too many requests...")
				sleep(30)
			else: return 0
	return 0

def fetch_summary(db, uid):  # Returns esummary results in the form of an XML string
	check_attempts()

	fetch_handle = Entrez.esummary(db=db, id=uid, retmode="xml")
	xml_doc = fetch_handle.read()
	fetch_handle.close()
	return xml_doc

def insert_gene(uid, xml_gene):  # The gene UID and name is inserted into our list of preliminary genes, assuming the gene was verified to have an association
	if xml_gene == 0: return  # Error occurred on fetch
	try:
		tree = ET.ElementTree(ET.fromstring(xml_gene))
		root = tree.getroot()
		symbol = (root.find(".//Gene-ref_locus")).text  # Location of gene symbol in the gene document
		preliminary_gene_uids[uid] = symbol
	except AttributeError:
		print("An unknown error occurred inserting gene... Invalid UID? - UID: " + uid)
		return

def verify_doc(db, uid, xml_doc):  # Checks to see that there is an association between the document and the search term
	verified_doc = False  # Flag to check if xml_doc (except gene documents) contains the search term
	tree = ET.ElementTree(ET.fromstring(xml_doc))
	root = tree.getroot()

	if db == "gene":
		for phenotype in root.findall(".//Gene-commentary/Gene-commentary_heading"):  # Location of heading items under the "Phenotypes" subsection in NCBI
			for term in mesh_terms:
				if term in phenotype.text.lower(): insert_gene(uid, xml_doc); return

	elif db == "medgen":
		for syn in root.findall(".//Name"):  # Location of search term synonyms
			for term in mesh_terms:
				if term in syn.text.lower():
					verified_doc = True
					break
			if verified_doc: break
		if verified_doc:
			for item in root.findall(".//Gene"):  # Loction of gene names
				if item.attrib['gene_id'] not in preliminary_gene_uids:  # Prevents adding the same gene
					parsed(item.attrib['gene_id'])
					xml_gene = fetch("gene", item.attrib['gene_id'])
					insert_gene(item.attrib['gene_id'], xml_gene)

	elif db == "gtr":
		with open("library/test_condition_gene.txt", "r") as in_file:  # This local file must be kept up to date
			for line in in_file:
				if uid in line:  # Find the UID in the document, it should always be there (if not, we need to update it from NCBI!!!)
					for term in mesh_terms:
						if term in line.lower():  # The line will contain condition information for the particular test, we match it with our MeSH terms
							verified_doc = True
							break
				if verified_doc: break

		if verified_doc:
			for item in root.findall(".//GeneID"):  # Location of gene UIDs
				if item.text not in preliminary_gene_uids:  # Prevents adding the same gene
					parsed(item.text)
					xml_gene = fetch("gene", item.text)
					insert_gene(item.text, xml_gene)

def seed_preliminary_gene_uids():  # Populates our list with a preliminary set of genes that have not been cross-checked
	if term == "": raise Exception("An empty search cannot be performed.")  # We must ensure that a term is initialized in R, not in Python

	# Perform a preliminary search in Gen
	print("Mining Gene...")
	gene_results = search("gene", (term+" "+organism))
	for uid in gene_results:
		if not parsed(uid):
			xml_gene = fetch("gene", uid)
			verify_doc("gene", uid, xml_gene)

	# Perform a preliminary search in MedGen
	print("Mining MedGen...")
	medgen_results = search("medgen", term)
	for uid in medgen_results:
		if not parsed(uid):
			xml_medgen = fetch_summary("medgen", uid)
			verify_doc("medgen", uid, xml_medgen)

	# Perform a preliminary search in GTR
	print("Mining GTR...")
	gtr_results = search("gtr", term)
	for uid in gtr_results:
		if not parsed(uid):
			xml_gtr = fetch_summary("gtr", uid)
			verify_doc("gtr", uid, xml_gtr)

def cross_check_gene(gene):  # The gene is crosschecked with PubMed to find literature with a co-occurence of the gene and term
	global term
	pubmed_results = search("pubmed", (term + " AND " + gene + "[Title/Abstract]"))  #####!!!!! BUG: Jun + Eng
	# Use of "AND" prevents dropping the gene term in the case that no article exists with the gene
	# Searching only the term would surely result in thousands of articles...
	# Example: "znf579 periodontitis"
	co_occurence_term, co_occurence_gene = False, False

	if len(pubmed_results) > 10: pubmed_uids = [pubmed_results[i] for i in range(10)]  # It should be found within 10 articles...
	else: pubmed_uids = pubmed_results.copy()

	print("Cross-checking: " + gene)
	for uid in pubmed_uids:
		xml_pubmed = fetch("pubmed", uid)  # PubMed articles can contain multiple genes, so do not add PubMed UIDs to parsed_uids
		tree = ET.ElementTree(ET.fromstring(xml_pubmed))
		root = tree.getroot()
		co_occurence_term, co_occurence_gene = False, False

		for abstract in root.findall(".//AbstractText"):
			if abstract.text != None:
				for term in mesh_terms:
					if term in abstract.text.lower():
						co_occurence_term = True
						break
				if gene.lower() in abstract.text.lower():
					co_occurence_gene = True
			if co_occurence_term == True and co_occurence_gene == True: break
		if co_occurence_term == True and co_occurence_gene == True: break

	if co_occurence_term == True and co_occurence_gene == True: return True
	else: return False

#---------- EXPORT DATA ----------#

def export_mesh_terms():
	with open("data/mesh_terms.txt", "w") as out_file:
		json.dump(mesh_terms, out_file)

def export_preliminary_gene_uids():
	with open("data/preliminary_gene_uids.json", "w") as out_file:
		json.dump(preliminary_gene_uids, out_file)	

def export_cross_checked_genes():
	with open("data/cross_checked_genes.csv", "w") as out_file:
		for gene in cross_checked_genes:
			out_file.write(gene + "\n")

def export_parsed_uids():
	with open("data/parsed_uids.json", "w") as out_file:
		json.dump(parsed_uids, out_file)		

#---------- DEBUG ----------#

#get_mesh_terms("alzheimer's disease")

#seed_preliminary_gene_uids()
#export_preliminary_gene_uids()

#for i in preliminary_gene_uids:
#	if crosscheck_gene(preliminary_gene_uids[i]):
#		crosschecked_gene_uids[i] = preliminary_gene_uids[i]
	
#export_crosschecked_gene_uids()