import streamlit as st
import time
from Primers.primer_functions import *


def page_primer():

    # ----------- Upload Fasta File --------
    fasta_file = st.file_uploader(label='Upload your Fasta file with genes',
                                  type='.fasta',
                                  help='Upload the file with the genes of interest')

    if fasta_file:
        fasta = load_fasta(fasta_file)

    # --- SideBar ---
    sidebar = st.sidebar

    if fasta_file:
        sidebar.write('Parameters')

        enzymes_to_use = list(sidebar.multiselect(label='Enzymes', options=[i for i in AllEnzymes]))
        sidebar.text('The order of selection matters')

        if enzymes_to_use:
            st.write('Order in vector:')
            draw_dna(enzymes_to_use)

        length_primers = sidebar.slider('Select the min/max length of primers',
                                        10, 35,
                                        (18, 24))

        temp_primers = sidebar.slider('Select the max difference of temperature between primers',
                                      0, 10,
                                      2)

        gc_primers = sidebar.slider('Select the min/max %GC for the primers',
                                    0, 100,
                                    (40, 60))

        end_gc = sidebar.checkbox('Allow only those ending with G or C', True)

        hairpin = sidebar.slider('Select the max value of Delta G (Gibbs free energy) to allow hairpins',
                                 -10, 10,
                                 -7)

        selfdimer = sidebar.slider('Select the max value of Delta G (Gibbs free energy) to allow Self Dimers',
                                   -10, 10,
                                   -7)

        crossdimer = sidebar.slider('Select the max value of Delta G (Gibbs free energy) to allow Cross Dimers',
                                    -10, 10,
                                    -7)

        pcr_limit = sidebar.slider('Select the minimum number of annealed nucleotides for in silico pcr',
                                   0, length_primers[1],
                                   13)

        run_app = sidebar.button('Run')

        enzymes_select = selected_enzymes(enzymes_to_use)
        # enzymes = selected_enzymes([AvrII, SpeI, AbsI, SgfI, BglII, AscI])
        enzymes_order = order_enzymes(enzymes_to_use)

    # --- Main ---
    if fasta_file and run_app:

        start = time.time()

        ids = fasta_ids(fasta)

        verify = verify_cuts(fasta, enzymes_select)

        primer = make_primers(fasta, length_primers[0], length_primers[1])

        pairs = make_pairs_primers(primer)

        unido = merge_tables(verify, pairs)

        pairs_enzymes = make_pairs_enzymes(unido, verify)

        order_filter = filter_by_order_vector_2(pairs_enzymes, enzymes_order)

        primers_with_enzymes = make_primers_enzymes_2(order_filter)

        temp = calculate_temp_2(primers_with_enzymes)

        filter_temp = filter_by_temp_2(temp, temp_max = int(temp_primers))

        gc_calc = calculate_GC_2(filter_temp)

        gc_filter = filter_by_GC_2(gc_calc, gc_primers[0], gc_primers[1])

        if end_gc:
            gc_filter = filter_by_end_GC_2(gc_filter)

        calc_hairpin = calculate_hairpin_2(gc_filter)

        filter_hairpin = filter_by_Hairpin_2(calc_hairpin, hairpin)

        self_dimer = calculate_selfdimer_2(filter_hairpin)

        filter_self_dimer = filter_by_Selfdimer_2(self_dimer, selfdimer)

        cross_calc = calculate_Crossdimer_2(filter_self_dimer)

        filter_cross = filter_by_Crossdimer_2(cross_calc, crossdimer)

        tamanho = calculate_len_2(filter_cross)

        pcr = pcr_in_silico(fasta, tamanho, pcr_limit)
        st.write(pcr)

        count = count_primers(pcr)

        count_gene = count_primers_for_each_gene(count, ids)

        count_all = count_global_primers(count, ids)
        st.write(count_all)

        group_pcr = group_table(pcr)

        excel = make_excel(pcr, group_pcr)

        st.download_button(
            label="Download Primers table",
            data=excel,
            file_name="Primers.xlsx",
            mime="application/vnd.ms-excel")

        end = time.time()
        total_time = end - start

        # st.write(f'A execução levou: {str(total_time)} segundos')
        # print(f'A execução levou: {str(total_time)} segundos')

        # -------------------------------------------------------------------
