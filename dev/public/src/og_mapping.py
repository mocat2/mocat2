#!/usr/bin/env python
"""
To see possible option run 'og_mapping.py -h'.

The aim of this program is to take environmental genes/proteins
("region") and figure out the orthologous groups (OG) they belong
to. For this the regions are blasted against STRING proteins. For each
STRING protein it is known which part of it contains sequences that
belong to a certain orthologous group. Depending on the overlap of
environmental and STRING sequences (and based on some heuristics) OGs
can be allocated to the environmental proteins.

A little drawing for your understanding:

--|==================================|-----------|===|--    Contigs with regions
   .                     .   .       .           .   .
   .                     .   .       .           .   .
   |=====================|   |================| |=========| STRING proteins
   .        .   .      .       .  .
   .        .   .      .       .  .
   <-------->   <------>       <-->                         OGs


                              String protein 1 - No COG
                             /
One environmental protein ---- Strint protein 2 - COG BLUB
                             \
                              String protein 3-+- COG XYZ
                                               |
                                               +- COG ABC
                                               |
                                               +- COG WHO
"""

import sys
from optparse import OptionParser

def main():
    """Run the mapping."""
    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(2)
    mapper = OGMapper()
    mapper.set_parameters()
    mapper.check_parameters()
    mapper.perform_mapping()
    mapper.generate_report()
        
class OGMapper(object):
    """See head of file for further description."""

    def __init__(self):
        self.count_factor = 1
        self.regions_and_mappings = {}
        self.unfiltered_og_mapping = []
        self.regions_best_hits = {} 
        self.regions_and_attributes = {}
        self.string_protein_og_mappings = {}        
        self.string_proteins_and_length = {}
        self.options = None

    def set_parameters(self):
        """Set command line paramters"""
        option_parser = OptionParser()
        option_parser.add_option(
           "-a", "--blast-algorithm", dest="blast_algorithm",
           help="algorithm of the BLAST program used. blastp or blastx " +
           "(default: blastp)", metavar="STRING", default="blastp")
        option_parser.add_option(
            "-f", "--blast-flavor", dest="blast_flavor",
            help="flavor of protein BLAST program used (default: NCBI)", 
            metavar="STRING", default="NCBI")
        option_parser.add_option(
            "-b", "--blast-file", dest="blast_file",
            help="input file with all BLAST hits", metavar="FILE")
        option_parser.add_option(
            "-t", "--best-blast-hit-file", dest="best_blast_hit_file",
            help="input file only with the best BLAST hits", metavar="FILE")
        option_parser.add_option(
            "-r", "--contig-region-mapping-file",
            dest="contig_region_mapping_file",
            help="contig region mapping file", metavar="FILE")
        option_parser.add_option(
            "-g", "--gene-orthgroup-file", dest="gene_orthgroup_file",
            help="STRING genes orthgroup file", metavar="FILE")
        option_parser.add_option(
            "-p", "--string-protein-properties-file", dest="string_poperties",
            help="file with string protein properties (e.g. length)", 
            metavar="FILE")
        option_parser.add_option(
            "-m", "--minimum-relevant-overlap-fraction",
            dest="minimum_relevant_overlap_fraction",
            help="COGs need to overlap by at least that much to become " +
            "incompatible (default = 0.5)", metavar="FLOAT", default=0.5)
        option_parser.add_option(
            "-s", "--minimum-bitscore", dest="min_bitscore",
            help="do not consider any COG-mappings below that bitscore " +
            "(default = 60)", metavar="FLOAT", default=60)
        (self.options, args) = option_parser.parse_args()

    def check_parameters(self):
        if (self.options.blast_algorithm == "blastx" and 
            self.options.contig_region_mapping_file == None):
            sys.stderr.write(
                'Error - You selected "blastx" as BLAST algorithm but you ' +
                'did not provide a contig-region-mapping file.\n')
            sys.exit(2)
        if (self.options.blast_algorithm == "blastp" and 
            not self.options.contig_region_mapping_file == None):
            sys.stderr.write(
                'You selected "blastp" as BLAST algorithm. ' +
                'In this case a contig-region-mapping file is not needed.\n')

    def perform_mapping(self):
        """Do the actual mapping."""
        self.__parse_genes_orthgroups_file()
        if self.options.blast_algorithm == "blastx": 
            self.__parse_contig_region_mapping_file()
        self.__parse_best_blast_hits_file()
        self.__parse_string_protein_length_file()
        self.__parse_blast_file_and_do_mapping()
        self.__filter_mappings()

    def __filter_mappings(self):
        """
        At this state we have a collection of region -> OG mappings
        based on the region -> STRING protein mappings and the STRING
        protein -> OG mappings. Now these mappings are filtered and
        evaluated.
        """
        for mapping in self.unfiltered_og_mapping:
            if self.__is_kog_or_twog(mapping): continue
            if mapping['bitscore'] < float(self.options.min_bitscore): continue
            if mapping['protein'] in self.regions_and_mappings:
                self._check_clash_with_prev_mappings(mapping)
            else:
                self.regions_and_mappings[mapping['protein']] = [mapping]

    def generate_report(self):
        """Print a report to standard output"""
        print("#Predicted_gene\tSTRING_gene\tCOG AlignmentStart_on_String_"
               "Protein\tAlignmentEnding_on_String_Protein\tSW_score\t"
               "Gene_length\tBEST:Best_STRING_Gene:Best_score")
        for region in sorted(self.regions_and_mappings.keys()):
            if region in self.regions_best_hits:
                best_hit_partner = self.regions_best_hits[region][0]
                best_hit_score = self.regions_best_hits[region][1]
            else:
                best_hit_partner = "NONE"
                best_hit_score = 0
            for mapping in self.regions_and_mappings[region]:
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tBEST:{7}:{8}".format(
                    region, mapping["STRING_protein"], mapping["OG"], 
                    mapping["start_STRING"], mapping["end_STRING"], 
                    mapping["bitscore"], self.__get_string_protein_length(
                        mapping["STRING_protein"]), best_hit_partner, 
                    best_hit_score))

    def __get_string_protein_length(self, protein):
        """Return the length of a given STRING protein"""
        return self.string_proteins_and_length[protein]

    def _check_clash_with_prev_mappings(self, mapping):
        """Test if there are collisions with previous mappings"""
        has_overlap = False
        # Important: Take a copy ([:]) of the returned list not the
        # list itself as the list might become short during the
        # process so not all items are checked.
        for prev_mapping in self._sort_prev_mappings_by_score(
            mapping['protein'])[:]: 
            # If there is an overlap with a previous sequence that is
            # large enough, the sequence might be discarded. All
            # previous mappings have to be checked even if one overlap
            # is found already.
            if (mapping['OG'] == prev_mapping['OG'] and
                self._has_relevant_overlap(mapping, prev_mapping)):
                has_overlap = True
                # If the overlap fraction is big enough to discard a
                # sequence, check if the previous sequence is
                # shorter. If this is the case this short sequence is
                # replaced by the new longer sequence. As a huge
                # sequence could have replaced more than one smaller
                # sequence - duplicated occurences have to be avoided
                # by searching before replacing.
                length_of_prev_mapping = self.__length_of_mapping(prev_mapping)
                length_of_cur_mapping = self.__length_of_mapping(mapping)
                if length_of_prev_mapping < length_of_cur_mapping:
                    # If the current big mapping is already in the
                    # list as it has already replaced another sequence,
                    # the previous sequence which is shorter is removed
                    # from the list.
                    if mapping in self.regions_and_mappings[
                        mapping['protein']]:
                        self.regions_and_mappings[mapping['protein']].remove(
                            prev_mapping)
                    # If the current big mapping is not yet in the
                    # list it replaces the previous smaller one.
                    else:
                        prev_mapping_index = self.regions_and_mappings[
                            mapping['protein']].index(prev_mapping)
                        self.regions_and_mappings[
                            mapping['protein']][prev_mapping_index] = mapping
                has_overlap = True
        # If there has been no overlap detected to any previous
        # mapping add the current mapping to the list
        if not has_overlap:
            self.regions_and_mappings[mapping['protein']].append(mapping)
        
    def _has_relevant_overlap(self, mapping, prev_mapping):
        """Test if the overlap to previous mapping is large enough"""
        _has_relevant_overlap = False
        overlap = float(self.__calc_overlap(
                mapping['start'], mapping['end'], 
                prev_mapping['start'], prev_mapping['end']))
        length_of_prev_mapping = self.__length_of_mapping(prev_mapping)
        length_of_cur_mapping = self.__length_of_mapping(mapping)
        if (overlap / length_of_prev_mapping >=
            float(self.options.minimum_relevant_overlap_fraction) or
            overlap / length_of_cur_mapping >=
            float(self.options.minimum_relevant_overlap_fraction)):
            _has_relevant_overlap = True
        return _has_relevant_overlap
    
    def __length_of_mapping(self, mapping):
        return float(mapping['end'] - mapping['start'])

    def _sort_prev_mappings_by_score(self, protein):
        """Sort the previous mappings by bitscore"""
        scores_and_mapping = []
        for mapping in self.regions_and_mappings[protein]:
            scores_and_mapping.append([mapping['bitscore'], mapping])
        scores_and_mapping.sort()
        scores_and_mapping.reverse()
        score_sorted_mapping = []
        for score, mapping in scores_and_mapping: 
            score_sorted_mapping.append(mapping)
        return score_sorted_mapping

    def __is_kog_or_twog(self, mapping):
        """Check if this is a KOG or a TWOG entry"""
        if mapping['OG'].startswith('KOG') or mapping['OG'].startswith('TWOG'):
            return True

    def __parse_contig_region_mapping_file(self):
        """
        Parse the contig to region mapping file. There are different
        line formats. One with 5 columns, one with 6 columns.
        """
        for split_line in self.__line_splitter(
            self.options.contig_region_mapping_file):
            contig, region, start, end, orientation = ['', '', '', '', '']
            if len(split_line) == 5: 
                contig, region, start, end, orientation = split_line
                origin = "default"
            elif len(split_line) == 6: 
                contig, region, start, end, orientation, origin = split_line
            self.regions_and_attributes[region] = {
                'start' : int(start), 'end' : int(end), 'origin' : origin,
                'orientation' : orientation}

    def __parse_string_protein_length_file(self):
        """Get the lenght of STRING proteins from the input file."""
        for split_line in self.__line_splitter(self.options.string_poperties):
            protein, length = split_line
            self.string_proteins_and_length[protein] = length

    def __parse_genes_orthgroups_file(self):
        """
        Extract the protein to orthgroup mapping and store them by
        STRING protein ids. One STRING protein can be mapped to
        multiple orthologous groups.
        """
        for split_line in self.__line_splitter(
            self.options.gene_orthgroup_file):
            og_id, protein, start, end = split_line[:4]
            if protein in self.string_protein_og_mappings:
                self.string_protein_og_mappings[protein].append(
                    {'OG' : og_id, 'start' : int(start), 'end' : int(end)})
            else:
                self.string_protein_og_mappings[protein] = [
                    {'OG' : og_id, 'start' : int(start), 'end' : int(end)}]

    def __parse_blast_file_and_do_mapping(self):
        """
        Parses the BLAST result file and performs the the OG mappings
        stepwise.
        """
        previous_protein = None
        cur_og_mappings = []
        # Got trought all region to STRING protein blast results
        for cur_hit_values in self.__get_blast_results(
            self.options.blast_file, self.options.blast_flavor):
            if cur_hit_values == None: continue
            cur_protein = cur_hit_values['region']
            cur_string_hit = cur_hit_values['STRING_protein']
            # If for a region the end of the list of all STRING hits is reached,
            # the collected OG mapping are processed
            if cur_protein != previous_protein and previous_protein != None:
                self.__process_prev_protein_results(
                    previous_protein, cur_og_mappings)
                cur_og_mappings = []
            previous_protein = cur_protein
            if not self.__hit_is_okay(cur_hit_values): continue
            # If a valid hit to a STRING protein is found, all valid
            # OG mappings of this STRING protein are added to the list
            # of OG mappings of the current region.
            self.__add_string_og_mappings_to_region(
                cur_hit_values, cur_og_mappings)
        self.__process_prev_protein_results(previous_protein, cur_og_mappings)

    def __add_string_og_mappings_to_region(
        self, cur_hit_values, cur_og_mappings):
        """ """
        for og_mapping in self.string_protein_og_mappings[
            cur_hit_values['STRING_protein']]:
            self.__add_mapping_to_region(
                cur_hit_values, og_mapping, cur_og_mappings)

    def __add_mapping_to_region(
        self, cur_hit_values, og_mapping, cur_og_mappings):
        """
        Compare the STRING-OG BLAST hit position with the STRING-environmental
        blast hit position.
        """
        og_id = og_mapping['OG']
        if not self.__overlap_good_enough(cur_hit_values, og_mapping): return
        mapping_start, mapping_end = self.__calc_mapping_positions(
            cur_hit_values, og_mapping)
        need_new_mapping_needed = self.__compare_with_previous_mappings(
            cur_hit_values, mapping_start, mapping_end, og_id, cur_og_mappings)
        if need_new_mapping_needed:
            cur_og_mappings.append(
                self.__get_mapping_values(
                    cur_hit_values, mapping_start, mapping_end, og_id))

    def __overlap_good_enough(self, cur_hit_values, og_mapping):
        """
        The overlap of the region with the STRING protein and the
        sequence of OG with the STRING protein have to overlap to
        certain level to make this a valid mapping. Return True is
        this is the case.
        """
        overlap = self.__calc_overlap(
            cur_hit_values['start_STRING'], cur_hit_values['end_STRING'],
            og_mapping['start'], og_mapping['end'])
        return overlap > 10 * self.count_factor

    def __get_mapping_values(
        self, cur_hit_values, mapping_start, mapping_end, og_id):
        """ """
        if self.options.blast_algorithm == "blastx": 
            strand = self.__get_strand(cur_hit_values['region'])
        else: strand = ''
        return {'bitscore' : cur_hit_values['bitscore'],
                'STRING_protein' : cur_hit_values['STRING_protein'],
                'start_STRING' : cur_hit_values['start_STRING'],
                'end_STRING' : cur_hit_values['end_STRING'],
                'start' : mapping_start,
                'end' : mapping_end,
                'OG' : og_id,
                'strand' : strand}

    def __calc_mapping_positions(self, cur_hit_values, og_mapping):
        """
        Calculate/estimate the position of the OG on the environmental
        protein. This is a rough heurisms - a BLAST or Smith/Waterson
        alignment would be more precise.
        """
        string_start = cur_hit_values['start_STRING']
        string_end = cur_hit_values['end_STRING']
        hit_start = cur_hit_values['start_region']
        hit_end = cur_hit_values['end_region']
        string_og_start = og_mapping['start']
        string_og_end = og_mapping['end']
        length_ratio = self.__calc_length_ratio(cur_hit_values)
        mapping_start = 0
        mapping_end = 0
        if self.options.blast_algorithm == "blastp":
            mapping_start = int(hit_start + (
                    float(string_og_start - string_start) / length_ratio))
            mapping_end = int(hit_end - (
                    float(string_end - string_og_end) / length_ratio))
        elif self.options.blast_algorithm == "blastx":
            strand = self.__get_strand(cur_hit_values['region'])
            if strand == "forward":
                mapping_start = int(hit_start + (
                        float(string_og_start - string_start) / length_ratio))
                mapping_end = int(hit_end - (
                        float(string_end - string_og_end) / length_ratio))
            else:
                mapping_start = int(hit_start + (
                        float(string_end - string_og_end) / length_ratio))
                mapping_end = int(hit_end - (
                        float(string_og_start - string_start) / length_ratio))
        if mapping_start < hit_start: mapping_start = hit_start
        if mapping_end > hit_end: mapping_end = hit_end
        return mapping_start, mapping_end

    def __process_prev_protein_results(self, protein, cur_og_mappings):
        """Once all OG mapping of a region are detected they are stored """
        for mapping_values in cur_og_mappings:
            self.unfiltered_og_mapping.append({
                'protein' : protein,
                'STRING_protein' : mapping_values['STRING_protein'],
                'OG' : mapping_values['OG'],
                'start' : mapping_values['start'],
                'end' : mapping_values['end'],
                'strand' : mapping_values['strand'],
                'start_STRING' : mapping_values['start_STRING'],
                'end_STRING' : mapping_values['end_STRING'],
                'bitscore' :mapping_values['bitscore']})

    def __compare_with_previous_mappings(
        self, cur_hit_values, mapping_start, mapping_end, og_id, 
        cur_og_mappings):
        """
        Compare the currently found OG mapping with the previous
        mapping on this region and see if this improved a previous
        one or lead to an new mapping.
        """
        for best_mapping_values in cur_og_mappings:
            overlap =  self.__calc_overlap(
                best_mapping_values['start'], best_mapping_values['end'],
                mapping_start, mapping_end)
            # Skip if this mapping does not seem to be a better
            # candidate for a previously found instance
            if (not og_id == best_mapping_values['OG'] or
                overlap < 20 * self.count_factor):
                continue
            if self.options.blast_algorithm == "blastx":
                strand = self.__get_strand(cur_hit_values['region'])
                if not strand == best_mapping_values['strand']: continue
            if mapping_start < best_mapping_values['start']:
                best_mapping_values['start'] = mapping_start
            if mapping_end > best_mapping_values['end']:
                best_mapping_values['end'] = mapping_end
            if cur_hit_values['bitscore'] > best_mapping_values['bitscore']:
                best_mapping_values['bitscore'] = cur_hit_values['bitscore']
                best_mapping_values['STRING_protein'] = cur_hit_values[
                    'STRING_protein']
            return False
        # There was no previous mapping that was similar to the
        # current one so start a new one.
        return True

    def __hit_is_okay(self, cur_hit_values):
        """ """
        if cur_hit_values['STRING_protein'] not in self.string_protein_og_mappings:
            sys.stderr.write('No OG mapping for STRING protein "%s"!\n' %
                             cur_hit_values['STRING_protein'])
            return False
        if cur_hit_values['end_STRING'] <= cur_hit_values['start_STRING']: 
            return False
        if cur_hit_values['start_region'] == cur_hit_values['end_region']: 
            return False
        return True

    def __get_blast_results(self, blast_result_file, blast_flavor):
        """Get the value from a BLAST line."""
        if blast_flavor == "NCBI": return self.__get_ncbi_blast_results(
            blast_result_file)
        if blast_flavor == "WU": return self.__get_wu_blast_results(
            blast_result_file)

    def __get_ncbi_blast_results(self, blast_result_file):
        """Extract NCBI BLAST values from a line"""
        for split_line in self.__line_splitter(blast_result_file):
            try:
                yield {'region' : split_line[0], 
                       'STRING_protein' : split_line[1],
                       'identity' : float(split_line[2]), 
                       'start_region' : int(split_line[6]), 
                       'end_region' : int(split_line[7]),
                       'start_STRING' : int(split_line[8]), 
                       'end_STRING' : int(split_line[9]),
                       'evalue' : float(split_line[10]), 
                       'bitscore' : float(split_line[11])}
            except IndexError:
                sys.stderr.write('Wrong NCBI BLAST output at line "%s"!' % 
                                 "\t".join(split_line))
                yield None

    def __get_wu_blast_results(self, blast_result_file):
        """Extract WU-BLAST values from a line"""
        for split_line in self.__line_splitter(blast_result_file):
            try:
                yield {'region' : split_line[0], 
                       'STRING_protein' : split_line[1],
                       'evalue' : float(split_line[2]), 
                       'bitscore' : float(split_line[4]),
                       'identity' : float(split_line[10]), 
                       'similarity' : float(split_line[11]),
                       'start_region' : int(split_line[17]), 
                       'end_region' : int(split_line[18]),
                       'start_STRING' : int(split_line[20]), 
                       'end_STRING' : int(split_line[21])}
            except IndexError:
                sys.stderr.write('Wrong WU-BLAST output at line "%s"!' % 
                                 "\t".join(split_line))
                yield None

    def __get_strand(self, protein):
        """Get the strand information"""
        return(self.regions_and_attributes[protein]['orientation'])

    def __calc_length_ratio(self, hit_values):
        """
        Calculate ration of the length of a String hit and a region
        hit
        """
        string_alignment_length = float(hit_values['end_STRING'] -
                                   hit_values['start_STRING'])
        region_alignment_lenth = float(hit_values['end_region'] -
                                  hit_values['start_region'])
        return(string_alignment_length / region_alignment_lenth)

    def __calc_overlap(self, start1, end1, start2, end2):
        """Calculate the overlap of two sequences"""
        overlap1 = end2 - start1
        overlap2 = end1 - start2
        # Todo. replace by min([overlap1,overlap2])
        overlap = overlap1
        if overlap2 < overlap1:
            overlap = overlap2
        return overlap

    def __line_splitter(self, file):
        """Return split lines and remove empty / commented lines."""
        for line in open(file):
            if line[0] in ['#', '\n']: continue
            split_line = line.split()
            yield split_line

    def __parse_best_blast_hits_file(self):
        """
        Parse the file containing the best blast hits of the
        regions against the STRING proteins
        """
        for split_line in self.__line_splitter(
            self.options.best_blast_hit_file):
            region = split_line[0]
            STRING_hit_protein = split_line[1]
            score = float(split_line[2])
            if score <= self.options.min_bitscore: continue
            self.regions_best_hits[region] = [STRING_hit_protein, score]

    # Todo: not used currently
    # def __calc_hit_positions(self, region_start, region_end, 
    #                          string_og_hit_start, string_og_hits_end, 
    #                          string_env_hit_start, string_env_hit_end, 
    #                          length_ratio, strand):
    #     """Calculate the postion of a hit"""
    #     hit_region_start = 0
    #     hit_region_end = 0
    #     if strand == "forward":
    #         hit_region_start = int(region_start + (float(
    #             string_og_hit_start - string_env_hit_start) / length_ratio))
    #         hit_region_end = int(region_end - (float(
    #             string_env_hit_end - string_og_hits_end) / length_ratio))
    #     else:
    #         hit_region_start = int(region_start + (float(
    #             string_env_hit_end - string_og_hits_end) / length_ratio))
    #         hit_region_end = int(region_end - (float(
    #             string_og_hit_start - string_env_hit_start) / length_ratio))
    #     if hit_region_start < region_start:
    #         hit_region_start = region_start
    #     if hit_region_end > region_end:
    #         hit_region_end = region_end
    #     return hit_region_start, hit_region_end
            
    # Todo not used?
    # def __sort_instances_by_position(self, instances, 
    #        instances_and_best_hits):
    #     instances_and_start_positions = []
    #     for instance in instances:
    #         instances_and_start_positions.append(
    #             [instances_and_best_hits[instance]['start'],
    #              instances_and_best_hits[instance]['end'], instance])
    #     instances_and_start_positions.sort()
    #     sorted_instances = []
    #     for start, end, instance in instances_and_start_positions:
    #         sorted_instances.append(instance)
    #     return sorted_instances

    # Old Version
    # def generate_report(self):
    #     """Print a report to standard output"""
    #     self.__parse_string_protein_length_file
    #     print ("#Predicted_gene\tSTRING_gene\tCOG AlignmentStart_on_String_"
    #            "Protein\tAlignmentEnding_on_String_Protein\tSW_score\t"
    #            "Gene_length")
    #     for region in sorted(self.regions_and_attributes.keys()):
    #         best_hit_partner = "NONE"
    #         best_hit_score = 0
    #         detection_type = "below_threshold"
    #         #attributes = self.regions_and_attributes[region]
    #         if self.regions_best_hits.has_key(region):
    #             detection_type = "homology"
    #             best_hit_partner = self.regions_best_hits[region][0]
    #             best_hit_score = self.regions_best_hits[region][1]
    #         # Discarded
    #         if not self.regions_and_mappings.has_key(region):
    #             continue
    #         for mapping in self.regions_and_mappings[region]:
    #             #start = mapping["start"] + attributes['start']
    #             #end = mapping["end"] + attributes['start']
    #             print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
    #                 region, mapping["STRING_protein"], mapping["OG"], 
    #                 mapping["start_STRING"], mapping["end_STRING"], 
    #                 mapping["bitscore"], self.__get_string_protein_length(
    #                     mapping["STRING_protein"]))

if __name__ == '__main__': main()
