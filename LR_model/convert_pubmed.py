# Author: Marcus Klang <marcus.klang@cs.lth.se>
#
# This script converts PubMed citation records into jsonl
#
# The script expects to be in the working directory with 
# all xml.gz files downloaded from:
# https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/

from lxml import etree
import gzip
import glob
from multiprocessing import Pool
import json
from tqdm import tqdm
import os

def iterate_xml(xmlfile):
    doc = etree.iterparse(xmlfile, events=('start', 'end'))
    _, root = next(doc)
    start_tag = None
    for event, element in doc:
        if event == 'start' and start_tag is None:
            start_tag = element.tag
        if event == 'end' and element.tag == start_tag:
            yield element
            start_tag = None
            root.clear()
            
def extract_data(node):
    def textof(path, default_val=None):
        el = node.find(path)
        if el is None:
            if default_val is not None:
                return default_val
            else:
                raise ValueError("No %s" % path)
        else:
            return el.text
        
    def abstract():
        elems = node.findall(".//AbstractText")
        abstracts = []
        for elem in elems:
            abstracts.append({"label": elem.attrib.get("Label", ""), "NlmCategory": elem.attrib.get("NlmCategory", ""), "Text": elem.text})
            
        return abstracts
        
    def mesh():
        mesh_root_el = node.find(".//MeshHeadingList")
        if mesh_root_el is None:
            return []
        
        meshes = []
        for mesh_heading in mesh_root_el.findall(".//MeshHeading"):
            mesh_heading_dict = {}
            
            descriptor = mesh_heading.find(".//DescriptorName")
            mesh_heading_dict["Name"] = descriptor.text
            mesh_heading_dict["UI"] = descriptor.attrib["UI"]
            mesh_heading_dict["MajorTopicYN"] =descriptor.attrib["MajorTopicYN"]
            
            qualifiers = []
            mesh_heading_dict["Qualifiers"] = qualifiers
            
            quals = mesh_heading.findall(".//QualifierName")
            for qual in quals:
                qual_dict = {}
                qual_dict["Name"] = qual.text
                qual_dict["UI"] = qual.attrib["UI"]
                qual_dict["MajorTopicYN"] = qual.attrib["MajorTopicYN"]
                qualifiers.append(qual_dict)
                
            meshes.append(mesh_heading_dict)
        
        return meshes
    
    def keywords():
        keyword_list = node.find(".//KeywordList[@Owner='NOTNLM']")
        if keyword_list is None:
            return []
        
        keywords = []
        for el in keyword_list.findall(".//Keyword"):
            keywords.append(el.text)
        
        return keywords
        
    
    PMID = textof(".//PMID")
    Title = textof(".//ArticleTitle")
    Abstract = abstract()
    Meshes = mesh()
    Keywords = keywords()
    PubDatePubMed = node.find(".//PubMedPubDate[@PubStatus='pubmed']")
    
    retval = {
        "PMID": PMID,
        "Title": Title,
        "Abstract": Abstract,
        "MeSH": Meshes,
        "Keywords": Keywords
    }
    
    if PubDatePubMed is not None:
        retval["PubDate"] = {
            "Year": int(PubDatePubMed.find(".//Year").text),
            "Month": int(PubDatePubMed.find(".//Month").text),
            "Day": int(PubDatePubMed.find(".//Day").text),
        }
    
    return retval

def process_file(path):
    target = os.path.basename(path)[0:-7] + ".jsonl.gz"
    with gzip.open(target, "wt") as fout:
        for node in iterate_xml(gzip.open(path, mode="rb")):
            json.dump(extract_data(node), fout)
            fout.write("\n")
    
    return path

if __name__ == "__main__":
    print("Starting conversion")
    with Pool(24) as pool:
        work = glob.glob("*.xml.gz")
        for res in tqdm(pool.imap_unordered(process_file, work), total=len(work)):
            pass
