import sys
import os
import re
from xml.dom import minidom
import pandas as pd
from Bio import Phylo
from bokeh import palettes
from ete3 import Tree
import json
import logging

class Metadata(object):
    """ This is the class for getting metadata for seqeunces.

    :param metadata_jsons: metadata of strains in JSON files
    :type metadata_jsons: list
    """

    def __init__(self, metadata_jsons):
        self.meta_df = self._load_metadata(metadata_jsons)

    def _load_metadata(self, metadata_files):
        mdf = pd.DataFrame()
        for file in metadata_files:
            if file.endswith(".json"):
                df = pd.read_json(file)
                df['taxonomy'] = df['Virus name'].str.replace('hCoV-19/', '')
                df['gender'] = df['Gender']
                df['age'] = df['Patient age']
                df['collection_date'] = df['Collection date']
                df['host'] = df['Host']
                df[['continent','country']] = df.Location.str.split(' / ',n=1, expand=True) 
                df['acc'] = df['Accession ID']
                df['sequencing_center'] = df['Submitting lab']
                df['coverage'] = df['Coverage']
            elif file.endswith(".tsv") or file.endswith(".txt"):
                df = pd.read_csv(file, sep="\t")
                df = df.append(mdf)

                # GISAID's metadata
                # strain  virus   gisaid_epi_isl  genbank_accession       date    region  country division        location        region_exposure country_exposure
                # division_exposure       segment length  host    age     sex     originating_lab submitting_lab  authors url     title   date_submitted
                if 'gisaid_epi_isl' in df:
                    df['taxonomy'] = df['strain']
                    df['acc'] = df['gisaid_epi_isl']
                    df['gender'] = df['sex']
                    df['collection_date'] = df['date']
                    df['continent'] = df['region']
                    df['submitting_lab'] = df['submitting_lab']
                    df['originating_lab'] = df['originating_lab']
                    df['authors'] = df['authors']
                    df['coverage'] = ""

                # GenBank's metadata
                # isolate accession col_date create_date Country: region
                if 'accession' in df:
                    import pycountry_convert as pc
                    
                    def country_to_continent(country_name):
                        country_continent_name = ""
                        try:
                            country_alpha2 = pc.country_name_to_country_alpha2(country_name)
                            country_continent_code = pc.country_alpha2_to_continent_code(country_alpha2)
                            country_continent_name = pc.convert_continent_code_to_continent_name(country_continent_code)
                        except:
                            pass

                        return country_continent_name
                    
                    df['acc'] = df['accession']
                    df['taxonomy'] = df['isolate'].str.replace('SARS-CoV-2/human/', '')
                    df['taxonomy'] = df['taxonomy'].str.replace('SARS-CoV-2/Homosapiens/', '')
                    df['taxonomy'] = df['taxonomy'].str.replace('SARS-CoV-2/', '')
                    df['collection_date'] = df['col_date']
                    df['country'] = df["Country: region"].str.split(r"\s?:\s?", expand=True)[0]
                    df['continent'] = df.apply(
                        lambda row: country_to_continent(row.country),
                        axis = 1
                    )
                    df['date_submitted'] = df['create_date']
                    df['date_submitted'] = pd.to_datetime(df['date_submitted'], errors='coerce')
                    df['date_submitted'] = df['date_submitted'].dt.strftime('%Y-%m-%d')
                    df['gender'] = ""
                    df['age'] = ""
                    df['originating_lab'] = ""
                    df['coverage'] = ""

            mdf = mdf.append(df)
        
        mdf = mdf.set_index('acc')
        return mdf

    def get_meta(self, acc):
        meta = {
            'taxonomy': "Unknown",
            'gender': "Unknown",
            'age': "Unknown",
            'collection_date': "Unknown",
            'host': "Unknown",
            'continent': "Unknown", 
            'country': "Unknown",
            'coverage': "Unknown",
            'submitting_lab': "Unknown",
            'originating_lab': "Unknown",
            'authors': "Unknown"
        }
        try:
            acc_series = self.meta_df.loc[[acc]].iloc[0]
            return acc_series[['taxonomy', 'gender', 'age', 'collection_date','host','continent','country','coverage','submitting_lab','originating_lab','authors']].fillna("None").to_dict()
        except:
            return meta

class PhyXML(object):
    """ This is the class for generating PhyD3 compatible tree file in extended phyloxml format.

    :param nwk_tree: The file location of the input newick tree
    :type nwk_tree: str
    :param xml_tree: The file location of the output XML file
    :type xml_tree: str
    :param metadata_jsons: metadata of strains in JSON files
    :type metadata_jsons: list
    """

    def __init__(self,
                 nwk_tree,
                 xml_tree,
                 assayseq,
                 Metadata_obj,
                 match_table_path,
                 dataSource):
        self.nwk_tree = nwk_tree
        self.xml_tree = xml_tree
        self.assayseq = assayseq
        self.acc_list = []
        self.acc_id_map = {}
        self.meta_json = {}
        self.meta = Metadata_obj
        self.dataSource = dataSource
        self.mt_df_orig = pd.read_csv(match_table_path, index_col=0)

        if self.dataSource == 'gisaid':
            cols = [c for c in self.mt_df_orig.columns if c and c.startswith('EPI_')]
            self.mt_df_orig = self.mt_df_orig[cols]
        elif self.dataSource == 'genbank':
            cols = [c for c in self.mt_df_orig.columns if c and not c.startswith('EPI_')]
            self.mt_df_orig = self.mt_df_orig[cols]
        
        self.mt_df = self._process_match_table()
        self.tree_obj = self._load_tree()
        self.dom = None
        self.colors = palettes.Category10[7]
        self.valid_result_df = pd.DataFrame()
        self.valid_result_orig_df = pd.DataFrame()
        self.taxa_colors = {}
        self.cont_colors = {
            'Asia': '0x1f77b4',
            'Europe': '0x2ca02c',
            'North America': '0xd62728', 
            'Africa': '0x8c564b',
            'Central America': '0xe377c2',
            'South America': '0x9467bd',
            'Oceania': '0x17becf',
            'Unknown': '0x7f7f7f',
        }
        self.stats = {
            "leaf_num": 0,
            "collapsed_genome_num": 0,
            "genome_num": 0,
            "genome_source": {
                "genbank": 0,
                "gisaid": 0
            },
            "collapsed_id": {},
            "nid_to_acc": {}
        }

    def _load_tree(self):
        # increase recursion limit
        sys.setrecursionlimit(100000)

        # proceeded tree
        proc_tree = self.nwk_tree
        return Tree(proc_tree)

    def purge_tree(self):
        """
        removing leaves not in overlap_table
        """
        # increase recursion limit
        sys.setrecursionlimit(100000)

        t = self.tree_obj
        df = self.mt_df_orig

        for node in t.traverse("postorder"):
            if node.name and ( node.name not in df or df[node.name][0]=='None'):
                node.delete()
        
        self.tree_obj = t


    def collapsing_nodes(self, 
                    collapse_node_by_branch=True,
                    collapse_node_by_pattern=True,
                    collapse_branch_len=0):

        def mean(array):
            return sum(array)/float(len(array))

        def cache_distances(tree):
            ''' precalculate distances of all nodes to the root''' 
            node2rootdist = {tree: 0}
            for node in tree.iter_descendants('preorder'):
                node2rootdist[node] = node.dist + node2rootdist[node.up]
                node.add_features(dist=node2rootdist[node])
            return node2rootdist


        def cache_heatmap(tree):
            ''' precalculate heatmap patterns of all nodes ''' 
            node2pattern = {tree: 'X'}
            df = self.mt_df

            for node in tree.iter_descendants('postorder'):
                if node.is_leaf():
                    ptn = 'X'
                    try:
                        name = node.name
                        if '=' in node.name:
                            name = node.name.split('=')[0]
                        text = "".join(df[name].astype(str).to_list())
                        if text:
                            ptn = text
                    except:
                        pass

                    node2pattern[node] = ptn
                    node.add_features(ptn = ptn)

                ptn = node2pattern[node]

                if node.up in node2pattern:
                    if node2pattern[node.up]:
                        # if the node has a different pattern, change concensus pattern to 'X'.
                        if node2pattern[node.up] != ptn:
                            node2pattern[node.up] = 'X'
                    else:
                        node2pattern[node.up] = ptn
                else:
                    # init pattern of internal node
                    node2pattern[node.up] = ptn

            return node2pattern


        def collapse_by_len(tree, min_dist):
            # cache the tip content of each node to reduce the number of times the tree is traversed
            node2tips = tree.get_cached_content()
            root_distance = cache_distances(tree)

            for node in tree.get_descendants('preorder'):
                if not node.is_leaf():
                    avg_distance_to_tips = mean([root_distance[tip]-root_distance[node]
                                                for tip in node2tips[node]])

                    if avg_distance_to_tips <= min_dist:
                        # rename
                        node.name += '='.join([tip.name for tip in node2tips[node]])
                        # label
                        node.add_features(collapsed_dist=True)


        def collapse_by_pattern(tree):
            # cache the tip content of each node to reduce the number of times the tree is traversed
            node2tips = tree.get_cached_content()
            node_pattern = cache_heatmap(tree)
            
            for node in tree.get_descendants('preorder'):
                if not node.is_leaf():
                    if node_pattern[node] != 'X':
                        # rename
                        node.name += '='.join([tip.name for tip in node2tips[node]])
                        # label
                        node.add_features(collapsed_ptn=True)


        # increase recursion limit
        sys.setrecursionlimit(100000)
        t = self.tree_obj
        #R = t.get_midpoint_outgroup()
        #t.set_outgroup(R)

        if collapse_node_by_branch:
            # label nodes that will be collapsed
            collapse_by_len(t, collapse_branch_len)
            # collapsed nodes are labeled, so you locate them and prune them
            for n in t.search_nodes(collapsed_dist=True):
                for ch in n.get_children():
                    ch.detach()

        if collapse_node_by_pattern:
            # label nodes that will be collapsed
            collapse_by_pattern(t)
            # collapsed nodes are labeled, so you locate them and prune them
            for n in t.search_nodes(collapsed_ptn=True):
                for ch in n.get_children():
                    ch.detach()

        # write to 
        t.write(outfile=f'{self.xml_tree}.temp_newick', format=5)
        proc_tree = f'{self.xml_tree}.temp_newick'            

        # converting to phyxml format
        Phylo.convert(proc_tree, 'newick', self.xml_tree, 'phyloxml')
        
        # parse phyloXML tree to DOM
        self.dom = minidom.parse(self.xml_tree)

    def _process_match_table(self):
        """
        Processing match table for tree pattern
        """
        # load validation tahble
        primer_df = self.mt_df_orig
        
        # transform data; 15 means assay failure
        def transform_func(x):
            if x == -1:
                return 3
            elif x != "None" and int(x) >= 3:
                return 3
            else:
                return x
        
        return primer_df.applymap(transform_func)


    def _create_text_node(self, tag, text, tag_attrs={}):
        x = self.dom.createElement(tag)
        for attr, val in tag_attrs.items():
            x.setAttribute(attr, val)
        txt = self.dom.createTextNode(text)
        x.appendChild(txt)
        return x

    def prepare_ext_phyloxml(self):
        itemlist = self.dom.getElementsByTagName('name')
        cc = list(self.colors)
        self.stats["leaf_num"] = len(itemlist)
        cnt = 0

        for elem in itemlist:
            acc = elem.firstChild.data
            # take care of reference genome
            if acc.startswith('NC_045512'):
                acc = 'NC_045512'

            # Collapsed node has a name that are concatenated with '_'
            collapsed_accs = ""
            acc = re.sub(r'(\d+)_', r'\1,', acc)

            if ',' in acc: #Collapsed node
                collapsed_accs = acc
                acc = collapsed_accs.split(',')[0]
                elem.firstChild.data = f'{acc}+'

            self.acc_list.append(acc)
            clade = elem.parentNode
            meta = self.meta.get_meta(acc)
            self.acc_id_map[acc] = str(cnt)
            self.stats['nid_to_acc'][str(cnt)] = acc
            cnt += 1

            taxa_name = ""
            if collapsed_accs:
                ids = collapsed_accs.split(',')
                country_cnt = {}
                continent_cnt = {}
                self.stats["collapsed_genome_num"] += len(ids)
                self.stats["genome_num"] += len(ids)
                self.stats["collapsed_id"][ids[0]] = ids

                for id in ids:
                    id_meta = self.meta.get_meta(id)
                    country = str(id_meta['country'])
                    continent = str(id_meta['continent'])
                    # counting the number of records for each country
                    if not country in country_cnt:
                        country_cnt[country] = 0
                    country_cnt[country] += 1
                    # counting the number of records for each continent
                    if not continent in continent_cnt:
                        continent_cnt[continent] = 0
                    continent_cnt[continent] += 1

                    if id.startswith("EPI"):
                        self.stats["genome_source"]["gisaid"] += 1
                    else:
                        self.stats["genome_source"]["genbank"] += 1

                # taxanomy name for collapsed node will be "country(count)" by default
                taxa_name = '; '.join([f"{c}({country_cnt[c]})" for c in country_cnt])

                if len(taxa_name)>40:
                    taxa_name = '; '.join([f"{c}({continent_cnt[c]})" for c in continent_cnt])
                if len(taxa_name)>40:
                    max_c = max(continent_cnt, key=lambda k: continent_cnt[k])
                    taxa_name = f"{max_c}({continent_cnt[max_c]}); total {len(ids)}"
                
                # if all countries belong to the same continent
                if len(continent_cnt) == 1:
                    meta['continent'] = list(continent_cnt.keys())[0]
                else:
                    meta['continent'] = 'Unknown'
            else:
                taxa_name = str(meta['taxonomy'])
                self.stats["genome_num"] += 1
                if acc.startswith("EPI"):
                    self.stats["genome_source"]["gisaid"] += 1
                else:
                    self.stats["genome_source"]["genbank"] += 1

            # add id
            x = self._create_text_node("id", self.acc_id_map[acc])
            clade.appendChild(x)

            # add taxonomy
            x = self.dom.createElement("taxonomy")
            c = self._create_text_node("code", taxa_name)
            x.appendChild(c)
            clade.appendChild(x)

            # add colortag likes <colortag>0xF52AD4</colortag>
            try:
                rgb = self.cont_colors[meta['continent']]
            except:
                rgb = self.cont_colors['Unknown']

            self.taxa_colors[acc] = {'taxa': taxa_name, 'rgb': rgb}

            # prepare JSON file
            idx = elem.firstChild.data
            self.meta_json[idx] = {}
            self.meta_json[idx].update(meta)
            self.meta_json[idx]["taxonomy"] = taxa_name

            # x = self._create_text_node("colortag", rgb)
            # clade.appendChild(x)
            # if collapsed_accs:
            #     self.meta_json[idx]["collapsed"] = collapsed_accs if len(collapsed_accs)<200 else f'{collapsed_accs[0:200]}...'
            #     x = self._create_text_node("property",
            #                                 collapsed_accs,
            #                                 {"ref": "Collapsed", "applies_to": "clade"})
            #     clade.appendChild(x)

            # for tag in meta:
            #     if tag != "taxonomy":
            #         val = str(meta[tag])
            #         x = self._create_text_node("property",
            #                                 str(meta[tag]),
            #                                 {"ref": tag, "applies_to": "clade"})
            #         clade.appendChild(x)

        # create graphs dom
        gs_dom = self.dom.createElement("graphs")
        root = self.dom.getElementsByTagName('phyloxml')[0]
        root.appendChild(gs_dom)

        # removing confidence elements
        elems = self.dom.getElementsByTagName('confidence')
        for elem in elems:
            parent = elem.parentNode
            parent.removeChild(elem) 

        return True

    def add_taxanomy_color(self):
        ts_dom = self.dom.createElement("taxonomies")
        for acc in self.acc_list:
            t_dom = self.dom.createElement("taxonomy")
            t_dom.setAttribute("code", self.taxa_colors[acc]['taxa'])

            n = self._create_text_node(
                "color", self.taxa_colors[acc]['rgb'])
            t_dom.appendChild(n)
            ts_dom.appendChild(t_dom)

        return ts_dom

    def add_graph_validation_heatmap(self, val_table, title):
        # load validation tahble
        primer_table = val_table
        primer_df = pd.read_csv(primer_table, index_col=0)
        self.valid_result_orig_df = self.valid_result_orig_df.append(primer_df)
        
        # transform data; 15 means assay failure
        def transform_func(x):
            if x == -1:
                return 3
            elif x != "None" and int(x) >= 3:
                return 3
            else:
                return x
        primer_df = primer_df.applymap(transform_func)

        # store to valid_result_df
        self.valid_result_df = self.valid_result_df.append(primer_df)
        df = primer_df

        # prepare heatmap
        g_dom = self.dom.createElement("graph")
        g_dom.setAttribute("type", "heatmap")
        # add graph title
        n = self._create_text_node("name", title)
        g_dom.appendChild(n)
        # create primer name in field nodes
        l = self.dom.createElement("legend")
        l.setAttribute("show", "1")

        for p_name in df.index.tolist():
            f = self.dom.createElement("field")
            n = self._create_text_node("name", p_name)
            f.appendChild(n)
            l.appendChild(f)

        # heatmap color settings
        g = self.dom.createElement("gradient")
        #n = self._create_text_node("name", "SpectralR")
        n = self._create_text_node("name", "phyHeatmap")
        c = self._create_text_node("classes", "4")
        g.appendChild(n)
        g.appendChild(c)
        l.appendChild(g)

        # add values
        d = self.dom.createElement("data")

        for acc in df.columns:
            if acc in self.acc_list:
                vs = self.dom.createElement("values")
                vs.setAttribute("for", self.acc_id_map[acc])
                for val in df[acc]:
                    v = self._create_text_node("value", str(val))
                    vs.appendChild(v)
                d.appendChild(vs)

        g_dom.appendChild(l)
        g_dom.appendChild(d)

        return g_dom

    def add_graph_bars(self, title):
        g_dom = self.dom.createElement("graph")
        g_dom.setAttribute("type", "multibar")

        n = self._create_text_node("name", title)
        g_dom.appendChild(n)

        # create a legend with fields
        l = self.dom.createElement("legend")
        l.setAttribute("show", "1")

        for text, color_code in zip(["Invalid (%)", "Valid (%)"], self.colors[1:3]):
            f = self.dom.createElement("field")
            n = self._create_text_node("name", text)
            f.appendChild(n)

            c = self._create_text_node("color", color_code)
            f.appendChild(c)
            l.appendChild(f)

        # add values
        d = self.dom.createElement("data")

        for acc in self.valid_result_df.columns:
            vs = self.dom.createElement("values")
            vs.setAttribute("for", self.acc_id_map[acc])
            tol_primer = len(self.valid_result_df[acc])

            primer_mismatch = pd.Series
            invalid_primer = 0

            try:
                primer_mismatch = self.valid_result_df[acc].value_counts()
            except:
                pass

            if 15 in primer_mismatch:
                invalid_primer += primer_mismatch[15]
            if "None" in primer_mismatch:
                invalid_primer += primer_mismatch["None"]

            # number of invalid primers
            v = self._create_text_node(
                "value", str(-invalid_primer/tol_primer*100))
            vs.appendChild(v)

            # number of valid primers
            v = self._create_text_node(
                "value", str((tol_primer-invalid_primer)/tol_primer*100))
            vs.appendChild(v)

            d.appendChild(vs)

        g_dom.appendChild(l)
        g_dom.appendChild(d)
        return g_dom

    def append_dom_to_tag(self, tag, dom):
        gs_dom = self.dom.getElementsByTagName(tag)[0]
        gs_dom.appendChild(dom)

    def writeXML(self):
        # minimizing XML
        fix = re.compile(r"\s{2,}")
        newXmlStr = re.sub(fix, '', self.dom.toxml())
        fix = re.compile(r'((?<=>)(\n[\t]*)(?=[^<\t]))|(?<=[^>\t])(\n[\t]*)(?=<)')
        newXmlStr = re.sub(fix, '', newXmlStr)
        # writing to output file
        with open(self.xml_tree, 'w') as file:
            file.write(newXmlStr)

    def get_stats(self):
        df1 = self.valid_result_orig_df.T.reset_index()
        df1['country'] = df1.apply(
            lambda d: self.meta.get_meta(d['index'])['country'] if not "/" in self.meta.get_meta(d['index'])['country'] else self.meta.get_meta(d['index'])['country'].split("/")[0].replace(" ","")
            , axis=1)
        df1['month'] = df1.apply(
            lambda d: self.meta.get_meta(d['index'])['collection_date'][:-3] if self.meta.get_meta(d['index'])['collection_date'].count('-')>1 else self.meta.get_meta(d['index'])['collection_date']
            , axis=1)
        df2 = df1.set_index(['country','month','index'])

        assay_df = pd.read_csv(
            self.assayseq, header=None, index_col=0, names=['assay','fp','rp','pb'], sep=" "
        )

        stats = {}
        def mm_stats(s):
            mmstats = {}
            for mm in s.index:
                if str(mm) != "None":
                    category = str(mm)
                    if int(mm) == 0:
                        continue
                    elif int(mm) < 0:
                        category = "8+"
                    elif int(mm) <= 3:
                        category = str(mm)
                    elif int(mm) <= 7:
                        category = "4-7"
                    else:
                        category = "8+"

                    if not category in mmstats:
                        mmstats[category] = int(s[mm])
                    else:
                        mmstats[category] += int(s[mm])
                else:
                    # category = "8+"
                    # if not category in mmstats:
                    #     mmstats[category] = int(s[mm])
                    # else:
                    #     mmstats[category] += int(s[mm])
                    pass
            return mmstats

        for assay in df2.columns:
            stats[assay] = {}
            stats[assay]['country'] = {}
            stats[assay]['month'] = {}
            stats[assay]['assay_sequence'] = {
                'forward_primer': assay_df.loc[assay,'fp'],
                'reverse_primer': assay_df.loc[assay,'rp'],
                'probe': assay_df.loc[assay,'pb']
            }
            
            for country in df2.index.levels[0]:
                s = df2.loc[[country],:][assay].value_counts()
                stat = mm_stats(s)
                if len(stat):
                    stats[assay]['country'][country] = stat

            idx = pd.IndexSlice
            for month in df2.index.levels[1]:
                s = df2.loc[idx[:,month],][assay].value_counts()
                stat = mm_stats(s)
                if len(stat):
                    stats[assay]['month'][month] = stat

        self.stats['assay_stats'] = stats

        return self.stats

    def writeMetaJson(self):
        # writing to output file
        with open(f'{self.xml_tree}.json', 'w') as file:
            json.dump(self.meta_json, file)

    def get_assay_list(self):
        return self.valid_result_df.index.tolist()

    def get_genome_list(self):
        return self.acc_list
