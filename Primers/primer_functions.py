from Bio import SeqIO
from Bio.Restriction import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from pydna.amplify import Anneal
from pydna.dseqrecord import Dseqrecord as Ds
from pydna.readers import read

import primer3

import pandas as pd
import numpy as np

import io

from io import StringIO

import streamlit as st

from dna_features_viewer import GraphicFeature, GraphicRecord

import matplotlib.pyplot as plt

import itertools


def load_fasta(fasta):
    stringio = StringIO(fasta.getvalue().decode("utf-8"))

    seqs = SeqIO.parse(stringio, 'fasta')

    return [seq for seq in seqs]


def fasta_ids(fasta):

  temp = []

  for gene in fasta:
    temp.append(gene.id)

  return temp


def selected_enzymes(list_enzymes):

  list_enzymes = RestrictionBatch(list_enzymes)
  return list_enzymes


def order_enzymes(list_enzymes):

  order_enzymes = dict()

  count = len(list_enzymes)

  for i in list_enzymes:

    order_enzymes[i] = count

    count -= 1

  return order_enzymes


def draw_dna(enzymes):
    count = 1
    features = []

    for i in enzymes:
        features.append(GraphicFeature(start=count,
                                       end=count + 6,
                                       strand=+1,
                                       color="#ffd700",
                                       linecolor='#ffd700',
                                       label=str(i)))
        count += 6

    record = GraphicRecord(sequence_length=count + 6, features=features)

    fig, ax = plt.subplots(figsize=(12, 2))

    record.plot(ax=ax)

    return st.pyplot(fig)


def verify_cuts(fasta, enzymes):

  prime_table = []

  for sequence in fasta:

    analyse = Analysis(enzymes, sequence.seq)

    for enzyme in analyse.without_site():

      prime_table.append([sequence.id, enzyme, enzyme.site])

  table_cuts = pd.DataFrame(prime_table, columns=['ID', 'Enzyme', 'Cut'])

  return table_cuts


def make_primers(fasta, range_min, range_max):
    record = []

    for seq in fasta:

        for lenght in range(range_min, range_max):
            fwd = str(seq.seq[:lenght])

            rev = str(seq.seq[-lenght:].reverse_complement())

            record.append({'ID': seq.id, 'Forward': fwd, 'Reverse': rev})

    table = pd.DataFrame(record)

    return table


def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)


def cartesian_product_multi(*dfs):
    idx = cartesian_product(*[np.ogrid[:len(df)] for df in dfs])
    return pd.DataFrame(
        np.column_stack([df.values[idx[:,i]] for i,df in enumerate(dfs)]))


def make_pairs_primers(tabela):

  temp = []

  for id, data in tabela.groupby('ID'):

    df1 = data[['ID', 'Forward']]
    df2 = data[['Reverse']]

    temp.append(cartesian_product_multi(*[df1, df2]))

  return pd.concat(temp).reset_index(drop=True).rename(columns = {0:'ID',	1: 'Forward', 2:'Reverse'})


def merge_tables(table_cuts, table_primers):

  return pd.merge(left = table_cuts, right = table_primers, on = 'ID')


def make_pairs_enzymes(tabela, verify):

  temp = []

  for id, data in tabela.groupby('ID'):

    df1 = verify[verify['ID'] == id]
    df2 = data[['Enzyme', 'Cut', 'Forward', 'Reverse']]

    temp.append(cartesian_product_multi(*[df1, df2]))

  return pd.concat(temp).reset_index(drop=True).rename(columns = {0: 'ID', 1: 'Enzyme_F', 2: 'Cut_F', 3: 'Enzyme_R', 4: 'Cut_R', 5: 'Forward', 6: 'Reverse'}).reset_index(drop = True)


def filter_by_order_vector_2(table, order):
    temp = []

    for id, data in table.groupby('ID'):
        temp.append(
            data[data.apply(lambda x: order[x['Enzyme_F']] > order[x['Enzyme_R']], axis=1)].reset_index(drop=True))

    return pd.concat(temp)


def make_primers_enzymes_2(table):

  temp = []

  for id, data in table.groupby('ID'):

    data['New_ID_F'] = data['Enzyme_F'].astype(str) + '_' + data['ID']

    data['New_ID_R'] = data['Enzyme_R'].astype(str) + '_' + data['ID']

    data['New_Forward'] = data['Cut_F'] + data['Forward']

    data['New_Reverse'] = data['Cut_R'].apply(lambda x: Seq(x).reverse_complement()).astype(str) + data['Reverse']

    temp.append(data)

  return pd.concat(temp)


def calculate_temp_2(table):

  temp = []

  for id, data in table.groupby('ID'):

    data['Temperature_F'] = data['New_Forward'].apply(lambda x: mt.Tm_NN(x))

    data['Temperature_R'] = data['New_Reverse'].apply(lambda x: mt.Tm_NN(x))

    data['Dif_Temp'] = abs(data['Temperature_F'] - data['Temperature_R'])

    temp.append(data)

  return pd.concat(temp)


def filter_by_temp_2(table, temp_max):

  temp = []

  for id, data in table.groupby('ID'):

    temp.append(data[(data['Dif_Temp'] <= temp_max)].reset_index(drop = True))

  return pd.concat(temp)


def calculate_GC_2(table):

  temp = []

  for id, data in table.groupby('ID'):

    data['GC_F'] = data['New_Forward'].apply(lambda x: gc_fraction(x)*100)
    data['GC_R'] = data['New_Reverse'].apply(lambda x: gc_fraction(x)*100)

    temp.append(data)

  return pd.concat(temp)


def filter_by_GC_2(table, percent_min, percent_max):
    temp = []

    for id, data in table.groupby('ID'):
        temp.append(data[((data['GC_F'] >= percent_min) & (data['GC_F'] <= percent_max)) &
                         ((data['GC_R'] >= percent_min) & (data['GC_R'] <= percent_max))].reset_index(drop=True))

    return pd.concat(temp)


def filter_by_end_GC_2(table):
    temp = []

    for id, data in table.groupby('ID'):
        temp.append(data[((data['New_Forward'].str[-1] == 'G') | (data['New_Forward'].str[-1] == 'C')) &
                         ((data['New_Reverse'].str[-1] == 'G') | (data['New_Reverse'].str[-1] == 'C'))].reset_index(
            drop=True))

    return pd.concat(temp)


def calculate_hairpin_2(table):

  temp = []

  for id, data in table.groupby('ID'):

    data['dG_Hairpin_F'] = data['New_Forward'].apply(lambda x: primer3.bindings.calc_hairpin(x).dg / 1000)
    data['dG_Hairpin_R'] = data['New_Reverse'].apply(lambda x: primer3.bindings.calc_hairpin(x).dg / 1000)

    temp.append(data)

  return pd.concat(temp)


def filter_by_Hairpin_2(table, dG):

  temp = []

  for id, data in table.groupby('ID'):

    temp.append(data[(data['dG_Hairpin_F'] >= dG) & (data['dG_Hairpin_R'] >= dG)].reset_index(drop = True))

  return pd.concat(temp)


def calculate_selfdimer_2(table):

  temp = []

  for id, data in table.groupby('ID'):

    data['dG_SelfDimer_F'] = data['New_Forward'].apply(lambda x: primer3.bindings.calc_homodimer(x).dg / 1000)
    data['dG_SelfDimer_R'] = data['New_Reverse'].apply(lambda x: primer3.bindings.calc_homodimer(x).dg / 1000)

    temp.append(data)

  return pd.concat(temp)


def filter_by_Selfdimer_2(table, dG):

  temp = []

  for id, data in table.groupby('ID'):

    temp.append(data[(data['dG_SelfDimer_F'] >= dG) & (data['dG_SelfDimer_R'] >= dG)].reset_index(drop = True))

  return pd.concat(temp)


def calculate_Crossdimer_2(table):

  temp = []

  for id, data in table.groupby('ID'):

    data['dG_CrossDimer'] = data.apply(lambda x: primer3.bindings.calc_heterodimer(x['New_Forward'], x['New_Reverse']).dg / 1000, axis = 1)

    temp.append(data)

  return pd.concat(temp)


def filter_by_Crossdimer_2(table, dG):

  temp = []

  for id, data in table.groupby('ID'):

    temp.append(data[data['dG_CrossDimer'] >= dG].reset_index(drop = True))

  return pd.concat(temp)


def calculate_len_2(table):
    temp = []

    for id, data in table.groupby('ID'):
        data['Lenght_F'] = data['New_Forward'].apply(lambda x: len(x))
        data['Lenght_R'] = data['New_Reverse'].apply(lambda x: len(x))

        temp.append(data)

    return pd.concat(temp)


def pcr_in_silico(fasta, table, limit):

  temp = []

  for j in table.iloc:

    for gene in fasta:

      if j['ID'] == gene.id:
        sequence = Ds(gene.seq)

        fwd = read(f'>{j.ID}\n{j.New_Forward}', ds=False)
        rev = read(f'>{j.ID}\n{j.New_Reverse}', ds=False)

        pcr = Anneal((fwd, rev), sequence, limit=limit)

        if len(pcr.products) == 1:

          temp.append(j)

  return pd.DataFrame(temp)


def count_primers(table):

  return dict(table.groupby('ID')['Enzyme_F'].count())


def count_primers_for_each_gene(group, fasta_ids):

  txt = []

  for fasta_id in fasta_ids:

    if fasta_id in group.keys():
      txt.append(f'The gene {fasta_id} has {group[fasta_id]} primers')

    else:
      txt.append(f'The gene {fasta_id} has no primers')

  return txt


def count_global_primers(group, fasta_ids):

  return (f'All genes: {len(fasta_ids)} | Genes with primers: {len(group)}')


def group_table(table):

  return table.groupby('ID')


def make_excel(table, group):

    buffer = io.BytesIO()

    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:

        for id, data in group:

            if id in table.ID.unique():
                data.to_excel(writer, sheet_name=f'{id}', index=False)

    #writer.save()

    return buffer