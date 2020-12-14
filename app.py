import json
import base64
import os
import dash
import pandas as pd
import dash_table
import subprocess
import dash_auth

import dash_html_components as html
import dash_core_components as dcc
import dash_uploader as du
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, title='AITumorBoard')
app.scripts.config.serve_locally = True

# Handle callback to component with id "fullband-switch"
app.config['suppress_callback_exceptions'] = True

text_style = {
    'color': "#506784",
    'font-family': "Open Sans"
}

VALID_USERNAME_PASSWORD_PAIRS = {
    'user001': '1234'
}

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')



def header_colors():
    return {
        'bg_color': '#0C4142',
        'font_color': 'white',
    }


app.layout = html.Div(id='alignment-body', className='app-body', children=[
    html.H1(" AITumorBoard"),
    html.H4(" AI-driven software platform for cancer comprehensive molecular analysis"),
    html.Div([
        html.Div(id='alignment-control-tabs', className='control-tabs', children=[
            dcc.Tabs(
                id='alignment-tabs', value='alignment-tab-select',
                children=[
                    dcc.Tab(
                        label='FASTQ data',
                        value='alignment-tab-select',
                        children=html.Div(className='control-tab', children=[
                            html.Div(className='app-controls-block', children=[

                                # Header "Status"
                                html.Div(id='container-button-basic3',
                                         className='fullwidth-app-controls-name',
                                         children='Status:'),

                                # Raw Status message:
                                html.Div(className='dashbio-loading', children=[
                                    html.Div(id='status-message',
                                             className='dashbio-loading',
                                             children='Provide the data and press "Calculate"',
                                             style={'width': 725, }),
                                ]),

                                # Raw Error message:
                                html.Div(className='dashbio-loading', children=[
                                    html.Div(id='err-message',
                                             className='dashbio-loading',
                                             children=' ',
                                             style={'width': 725, }),
                                ]),

                                html.Div(
                                    className='fullwidth-app-controls-name',
                                    children="Select the Patient disease:"
                                ),
                                dcc.Dropdown(
                                    id='alignment-dropdown',
                                    options=[
                                        {
                                            'label': 'Breast cancer',
                                            'value': 'dataset1'
                                        },
                                        {
                                            'label': 'Glioma',
                                            'value': 'dataset2'
                                        },
                                        {
                                            'label': 'Colorectal cancer',
                                            'value': 'dataset3'
                                        },
                                    ],
                                    value='dataset3',
                                )
                            ]),

                            # Provide links
                            html.Div(className='app-controls-block', children=[
                                html.Div(className='fullwidth-app-controls-name',
                                         children="Put links for the Patient DNA sequence data:"),

                                # window for patient name:
                                html.Div(dcc.Input(id='input-seq-1', type='text')),
                                html.Div(dcc.Input(id='input-seq-2', type='text')),

                                html.Div(id='invitation-2',
                                         className='fullwidth-app-controls-name',
                                         children='Or upload the data from local storage:'),


                                du.configure_upload(app, "/mnt/volume_nyc3_01/input", use_upload_id=False, upload_api=None),

                                html.Div([
                                    du.Upload(
                                        id='alignment-file-upload',
                                        filetypes=['csv', 'zip', 'gz', 'fq', 'fastq'],
                                        #max_files=2,
                                        text='Drag and drop FASTQ files here to upload!',
                                        text_completed='Uploaded: ',
                                        cancel_button=True,
                                        pause_button=True,
                                        max_file_size=25000,
                                        default_style=None,
                                    ),
                                ])


                            ]),

                            html.Div(id='container-button-basic',
                                     className='fullwidth-app-controls-name',
                                     children='Specify the Patient name:'),

                            # window for patient name:
                            html.Div(dcc.Input(id='input-name', type='text')),

                            # Start calc
                            html.Div(className='app-controls-block', children=[
                                html.Div(id='container-button-download2',
                                         className='fullwidth-app-controls-name',
                                         children='Start calculation if ready:'),

                                html.Button('Calculate', className='control-calculate', id='calculate', n_clicks=0),
                            ]),

                            # Interval module refreshing the app:
                            dcc.Interval(id='interval', interval=1 * 10000, n_intervals=0),

                        ])
                    ),
                    dcc.Tab(
                        label='Calculated',
                        value='control-tab-select2',
                        children=html.Div(className='control-tab', children=[

                            html.Div(
                                className='fullwidth-app-controls-name',
                                children="Select the calculated Patient:"
                            ),
                            dcc.Dropdown(
                                id='calculated-dropdown',
                                options=[

                                {'label' : na.split('.csv')[0], 'value': na} for na in os.listdir('/mnt/volume_nyc3_01/done/')

                                ],
                            )
                        ]),
                    ),
                    dcc.Tab(
                        label='About',
                        value='what-is',
                        children=html.Div(className='control-tab', children=[
                            html.H4(
                                className='what-is',
                                children='What is AI-TumorBoard?'
                            ),
                            html.P(
                                """
                                The AI-TumorBoard is a tool for comprehensive tumor 
                                genomic profiling. It allows to process genomic 
                                sequences from a FASTQ file to get actionable events
                                for targeted and immune therapies, and see a full
                                molecular profile of the tumor, to check which
                                pathways are damaged in cancer cells and in what
                                degree. It also provides populational statistics 
                                about similar cases from TCGA database.
                                It serves to oncology clinicians, cancer research 
                                scientists, and clinical research professionals.
                                To get the result you need to fill the fields on
                                the Data tab, and press Get result. Please scroll 
                                down to see the results.
                                After the calculation you can print the report.

                                """
                            ),
                            html.P(
                                """
                                Note that the AlignmentChart only returns a chart of
                                the sequence, while AlignmentViewer has integrated
                                controls for colorscale, heatmaps, and subplots allowing
                                you to interactively control your sequences.
                                """
                            ),
                            html.P(
                                """
                                Read more about the company here:
                                http://bioalg.com/
                                Contact us at:
                                info@bioalg.com
                                """
                            ),
                        ])
                    ),

                ],
            ),
        ]),
    ]),

    # table

    # status message:
    html.Div(className='dashbio-loading', children=[
        html.Div(id='result-field',
                 className='dashbio-loading',
                 children=' ',
    style={'width': 925,}),
    ]),


    html.Div(id='result-print'),

    # Result field
    html.Div(id='result-print2',
        style={'width': 725,}),

    # Result field as a pre-loaded table:



    html.Div(id='filename-print'),
    dcc.Store(id='alignment-data-store'),

])

auth = dash_auth.BasicAuth(
    app,
    VALID_USERNAME_PASSWORD_PAIRS
)

# Functions

#UPLOAD_DIRECTORY = "/data/input"


# Large files - special library


def action_at_click(se1, se2, value):
    """Calculation of full cycle: FASTQ processing"""

    # Put in file info that we are busy:
    semaphore.lock()

    # Firstly clear protocol file nohup:
    command00 = '> nohup.out'
    process = subprocess.Popen([command00], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    process.wait()

    # Read input data (3 items)
    pat = str(value)

    input_fastq1 = pat + '_raw_fastq_1.fq.gz'
    input_fastq2 = pat + '_raw_fastq_2.fq.gz'

    # Downloading:

    # Choose the way depending whether the link contains google:
    # **

    if 'google' in str(se1) or 'google' in str(se2):
        # for Google drive downloading:
        wget_first_part = 'wget --no-check-certificate -r "https://drive.google.com/uc?export=download&id='
        input_link1 = str(se1).split('https://drive.google.com/file/d/')[1].split('/view?usp=sharing')[0]
        input_link2 = str(se2).split('https://drive.google.com/file/d/')[1].split('/view?usp=sharing')[0]
        # with regex parse the string to get ID:
        # **
    else:
        # for AWS or Yandex downloading:
        wget_first_part = 'wget "'
        input_link1 = str(se1)
        input_link2 = str(se2)

    wget_second_part = '" -P /mnt/volume_nyc3_01/input/ -O /mnt/volume_nyc3_01/input/'

    command_download1 = wget_first_part + input_link1 + wget_second_part + input_fastq1
    command_download2 = wget_first_part + input_link2 + wget_second_part + input_fastq2

    # If we need to block the download
    process01 = subprocess.Popen([command_download1], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    process02 = subprocess.Popen([command_download2], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    process01.wait()
    process02.wait()

    # Check the size of input folder:
    # example:
    # sum(os.path.getsize(f) for f in os.listdir('.') if os.path.isfile(f))
    # if more 1 GB - go further

    path = '/mnt/volume_nyc3_01/input/'

    inputsize = sum(
        os.path.getsize(path + f) for f in os.listdir('/mnt/volume_nyc3_01/input/') if os.path.isfile(path + f))

    if inputsize < 1000000000:
        semaphore.unlock()
    else:

        # Start from pythonscipt:

        # fast to sam
        # command1 = 'bwa/bwa mem -t 8 /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta' + ' /mnt/volume_nyc3_01/input/' + input_fastq1 + ' /mnt/volume_nyc3_01/input/' + input_fastq2 + ' > /mnt/volume_nyc3_01/output/' + pat + '.sam'

        # if we not have sorted reads:
        # command1 = 'bowtie2 --very-sensitive --no-discordant --no-mixed -X 2000 -x /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta -1 /mnt/volume_nyc3_01/input/' + input_fastq1 + '-2 /mnt/volume_nyc3_01/input/' + input_fastq2 + ' -p 4 -S /mnt/volume_nyc3_01/output/' + pat + '.sam'
        # process = subprocess.Popen([command1], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        # process.communicate()
        # process.wait()

        # pipe - sample:
        # bwa mem -t 8 genome.fa reads.fastq | samtools sort -@8 -o output.bam -

        # this command not need because of pipe:
        # sam to bam
        # command2 = 'samtools view -@ 8 -S /mnt/volume_nyc3_01/output/' + pat + '.sam -b -o /mnt/volume_nyc3_01/output/' + pat + '.none.bam'
        # process2 = subprocess.Popen([command2], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        # process2.wait()

        # combine with alignment, since we don't need sam 20+GB
        # sort 1 time
        # command3 = 'samtools sort -@ 8 /mnt/volume_nyc3_01/output/' + pat + '.none.bam -o /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.bam'
        command3 = 'bwa/bwa mem -t 8 /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta /mnt/volume_nyc3_01/input/' + input_fastq1 + ' /mnt/volume_nyc3_01/input/' + input_fastq2 + ' | samtools sort -@ 8 -o /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.bam'
        process = subprocess.Popen([command3], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # fusion detection
        command35 = 'genefuse/genefuse -r /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta -f genefuse/genes/druggable.hg38.csv -1 /mnt/volume_nyc3_01/input/' + input_fastq1 + ' -2 /mnt/volume_nyc3_01/input/' + input_fastq2 + ' -h /mnt/volume_nyc3_01/output/GeneFuse_' + pat + '_report.html -t 8 > /mnt/volume_nyc3_01/output/GeneFuse_' + pat + '_result'
        process = subprocess.Popen([command35], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        # not wait for finishing this process, - let new be in parallel

        # read groups
        command4 = 'java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups -I /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.bam -O /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.bam -SO coordinate -RGID id -RGLB library -RGPL ILLUMINA -RGPU machine -RGSM OncoUnite'
        process = subprocess.Popen([command4], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # delete 01
        command45del = 'rm /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.bam'
        process = subprocess.Popen([command45del], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # sort bam
        command5 = 'java -jar picard/build/libs/picard.jar ReorderSam -I /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.bam -O /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.bam -SD /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta'
        process = subprocess.Popen([command5], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # delete 02
        command55del = 'rm /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.bam'
        process = subprocess.Popen([command55del], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # no duplicates
        command6 = 'java -jar picard/build/libs/picard.jar MarkDuplicates -I /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.bam -O /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.dedupped.bam -M output.metrics -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT'
        process = subprocess.Popen([command6], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # delete 03
        command65del = 'rm /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.bam'
        process = subprocess.Popen([command65del], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # quality
        command7 = 'gatk-4.1.8.1/gatk SplitNCigarReads -I /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.dedupped.bam -O /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.dedupped.split.bam -R /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta -skip-mq-transform true'
        process = subprocess.Popen([command7], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # Get VCF with Mutect:
        command8 = 'gatk-4.1.8.1/gatk --java-options "-Xmx12G" Mutect2 -R /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta -I /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.dedupped.split.bam -O /mnt/volume_nyc3_01/output/' + pat + '_mutect2_unfiltered.vcf'
        # --native-pair-hmm-threads 8'
        process = subprocess.Popen([command8], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # MSI detection using msisensor:
        command85 = 'msisensor2/msisensor2 msi -M msisensor2/models_hg38 -t /mnt/volume_nyc3_01/output/' + pat + '.none.sorted.RG.RS.dedupped.split.bam -o /mnt/volume_nyc3_01/output/msisensor2_output_' + pat + '.tumor.prefix'
        process = subprocess.Popen([command85], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        # not wait also

        # Filter VCF with Mutect2 function:
        command9 = 'gatk-4.1.8.1/gatk FilterMutectCalls -R /mnt/volume_nyc3_01/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta -V /mnt/volume_nyc3_01/output/' + pat + '_mutect2_unfiltered.vcf -O /mnt/volume_nyc3_01/output/' + pat + '_mutect2_filtered.vcf'
        process = subprocess.Popen([command9], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # Additional: filter for germline possibility < 50%:
        command10 = 'bcftools view -i "GERMQ<50" /mnt/volume_nyc3_01/output/' + pat + '_mutect2_filtered.vcf > /mnt/volume_nyc3_01/output/' + pat + '_mutect2_filtered_somatic.vcf'
        process = subprocess.Popen([command10], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # VEP annotate for AA change: (on NOT somatic and take full)
        command11 = 'ensembl-vep/vep -i /mnt/volume_nyc3_01/output/' + pat + '_mutect2_filtered_somatic.vcf -o /mnt/volume_nyc3_01/output/vep_output_' + pat + '_mutect2_filtered.txt -fork 8 -offline'
        process = subprocess.Popen([command11], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        process.wait()

        # end pythonscript

        # Interp function call for VEP
        interp2(pat)

        # Then clear the output folder:
        # rm -r /mnt/volume_nyc3_01/output/*

        command98_output_clear = 'rm -r /mnt/volume_nyc3_01/output/*'
        process = subprocess.Popen([command98_output_clear], shell=True, stdin=None, stdout=None, stderr=None,
                                   close_fds=True)
        process.wait()

        command99_input_clear = 'rm -r /mnt/volume_nyc3_01/input/*'
        process = subprocess.Popen([command99_input_clear], shell=True, stdin=None, stdout=None, stderr=None,
                                   close_fds=True)
        process.wait()

        semaphore.unlock()


def interp2(pat):
    """Interp 2 - for only VCF-vep"""
    # Start interp2

    # 1. Select VEP result to interpret:

    veps = pd.read_csv('/mnt/volume_nyc3_01/output/vep_output_' + pat + '_mutect2_filtered.txt', sep='\t', skiprows=40)

    veps = veps[veps['Protein_position'] != '-']

    veps = veps[veps['Amino_acids'].str.contains('/')]

    veps = veps[['Gene', 'Protein_position', 'Amino_acids']]

    # Split the column with AA

    veps = veps.assign(AA1=veps['Amino_acids'].str.split('/').str[0])
    veps = veps.assign(AA2=veps['Amino_acids'].str.split('/').str[1])

    # Two columns for connecting with oncoKB:

    veps['Alterations'] = veps['AA1'] + veps['Protein_position'] + veps['AA2']
    veps['Alter_f_Match'] = veps['AA1'] + veps['Protein_position']

    del veps['Protein_position']
    del veps['Amino_acids']
    del veps['AA1']
    del veps['AA2']

    # 2. Get Hugo names. Open the converter

    ensg = pd.read_csv('db/ENSG_converter.txt', sep='\t')

    ensg = ensg[['Approved symbol', 'Approved name', 'Chromosome', 'Ensembl ID(supplied by Ensembl)']]

    ensg.dropna(subset=['Ensembl ID(supplied by Ensembl)'], inplace=True)

    ensg.rename(columns={'Ensembl ID(supplied by Ensembl)': 'Gene'}, inplace=True)
    ensg.rename(columns={'Approved name': 'Short_Annotation'}, inplace=True)

    veps_hugo = veps.merge(ensg, on=['Gene'], how='inner')  # inner left

    del veps_hugo['Gene']
    veps_hugo.rename(columns={'Approved symbol': 'Gene'}, inplace=True)

    # 2.5 Add mutations effect from Cosmic:
    cos = pd.read_csv('db/cosmic_interp.csv', index_col=0)
    cos.rename(columns={'Gene name': 'Gene'}, inplace=True)

    # 1. First step: directly on full AA match:

    # to remove to COSMIC db preparation:
    cos['Alterations'] = cos['Alterations'].str.replace('=', '*')
    cos['Alterations'] = cos['Alterations'].str.replace('-', '*')

    veps_hugo['Alterations'] = veps_hugo['Alterations'].str.replace('=', '*')
    veps_hugo['Alterations'] = veps_hugo['Alterations'].str.replace('-', '*')

    veps_hugo_cos_direct_aa = veps_hugo.merge(cos, on=['Gene', 'Alterations'], how='inner')  # inner left
    del veps_hugo_cos_direct_aa['Alter_f_Match_x']
    del veps_hugo_cos_direct_aa['Alter_f_Match_y']

    # veps_hugo_cos_direct_aa = veps_hugo_cos_direct_aa[veps_hugo_cos_direct_aa['FATHMM prediction'] == 'PATHOGENIC']
    veps_hugo_cos_direct_aa.drop_duplicates(subset=['Gene', 'Alterations'], inplace=True)

    # Set trust = 1:
    veps_hugo_cos_direct_aa['Trust'] = 1

    # Rename and del columns before concat:
    veps_hugo_cos_direct_aa.rename(columns={'Alterations': 'Alteration'}, inplace=True)

    # 2. Second step: On match without last aminoacid *:

    veps_hugo['User_Alt'] = veps_hugo['Alterations']
    cos['Drug_Indicator'] = cos['Alterations']

    veps_hugo_cos_wth_last = veps_hugo.merge(cos, on=['Gene', 'Alter_f_Match'], how='inner')  # inner left
    # veps_hugo_cos_wth_last = veps_hugo_cos_wth_last[veps_hugo_cos_wth_last['FATHMM prediction'] == 'PATHOGENIC']
    veps_hugo_cos_wth_last.drop_duplicates(subset=['Gene', 'Alter_f_Match'], inplace=True)
    del veps_hugo_cos_wth_last['Alterations_x']
    del veps_hugo_cos_wth_last['Alterations_y']

    # Minus those who were in previous (full) match:

    veps_hugo_2nd_grade = veps_hugo_cos_wth_last[
        ~veps_hugo_cos_wth_last['User_Alt'].isin(veps_hugo_cos_direct_aa['Alteration'])]

    # Add trustness:

    veps_hugo_2nd_grade['Trust'] = 2
    veps_hugo_2nd_grade.loc[veps_hugo_2nd_grade['Drug_Indicator'].str.contains('\*'), 'Trust'] = 1

    # Rename and del cols before concat:
    del veps_hugo_2nd_grade['Alter_f_Match']
    del veps_hugo_2nd_grade['Drug_Indicator']
    veps_hugo_2nd_grade.rename(columns={'User_Alt': 'Alteration'}, inplace=True)

    # Concat 2 trust by Cosmic categories:
    cosmic_annotated = pd.concat([veps_hugo_cos_direct_aa, veps_hugo_2nd_grade], sort=True)

    # Sort
    cosmic_annotated.sort_values(by=['Trust', 'Gene'], inplace=True)

    # Left only with Trust = 1 (playable):
    cosmic_annotated = cosmic_annotated[cosmic_annotated['Trust'] == 1]

    # Temporry delete (since all  == 1):
    del cosmic_annotated['Trust']

    # Play with Fathmm score:

    # 3. Merge with drugs:

    act = pd.read_csv('db/oncokb_biomarker_drug_associations-2.tsv', sep='\t')
    ''' 'Level', 'Gene', 'Alterations', 'Tumor Type', 'Drugs' '''

    # Slice only indicators with Oncogenic mutation (before excluding):
    act_oncogenic = act[act['Alterations'].str.contains('Mutations')]  # 95 positions

    # Explode multiple alterations in string:
    act = act.assign(Alterations=act['Alterations'].str.split(', ')).explode('Alterations')

    # Add * to alterations without last letter:
    act.loc[act['Alterations'].str[-1].isin(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']), 'Alterations'] += '*'

    # Exclude
    search_list = '|'.join(
        ['del', 'ins', 'Fusion', 'dup', 'Amplification', 'Instability', 'Exon', 'Duplication', 'Wildtype', 'Burden',
         'Mutations', 'Deletion'])
    act = act[~act['Alterations'].str.contains(search_list)]

    # Two ways of drugs alterations:

    # 1. Directly on full gene and alteration:

    act['Alteration'] = act['Alterations']
    act_on_gene_alt = cosmic_annotated.merge(act, on=['Gene', 'Alteration'], how='inner')

    # 2. If AA with '*' symbol - match on Alter_f_match:

    act_no_last_aa = act[act['Alterations'].str.contains('\*')]  # 39

    # Add for soft search (with *):
    act_no_last_aa['Alter_f_Match'] = act_no_last_aa['Alterations'].str[:-1]
    del act_no_last_aa['Alteration']
    cosmic_annotated['Alter_f_Match'] = cosmic_annotated['Alteration'].str[:-1]
    act_on_gene_no_last = cosmic_annotated.merge(act_no_last_aa, on=['Gene', 'Alter_f_Match'], how='inner')
    del act_on_gene_no_last['Alter_f_Match']
    del cosmic_annotated['Alter_f_Match']

    # 3. If Oncogenic mutation - then connect with any pathogenic mutation:

    # Left only pathogenic:
    cosmic_annotated_patho = cosmic_annotated[cosmic_annotated['FATHMM prediction'] == 'PATHOGENIC']

    act_on_gene_oncogenic = cosmic_annotated_patho.merge(act_oncogenic, on=['Gene'], how='inner')  # 12 rows

    # Concat all 3 results:

    act_res = pd.concat([act_on_gene_alt, act_on_gene_no_last, act_on_gene_oncogenic], sort=True)

    # Rename user alteration
    # act_on_gene_oncogenic.rename(columns={'Alteration': 'Molecular_Event'}, inplace=True)
    # act_on_gene_oncogenic.rename(columns={'Alterations': 'Drug_Indicator'}, inplace=True)

    # Modify: group by

    act_on_gene_oncogenic_gr_drug = act_res.groupby(['Gene', 'Tumor Type', 'Alteration'])['Drugs'].apply(
        ', '.join).reset_index()  # .groupby('Gene')['Drugs'].apply(', '.join).reset_index()

    # del act_on_gene_oncogenic_gr_event['Drugs']

    def left_unique_in_cell(gene_diseases_drug, input_col, output_col):
        """ delete duplicates in cell (list) """
        gene_diseases_drug['temp'] = gene_diseases_drug[input_col].str.split(',')

        gene_diseases_drug[output_col] = ''

        for ind in gene_diseases_drug.index:
            # get the num of items:
            list_drugs = gene_diseases_drug.at[ind, 'temp']
            list_drugs = [x.strip() for x in list_drugs]
            gene_diseases_drug.at[ind, output_col] = list(set(list_drugs))

        del gene_diseases_drug['temp']
        del gene_diseases_drug[input_col]

    def take_off_brackets_quotes(df, col):
        df[col] = df[col].astype(str)
        df[col] = df[col].str.replace('[', '')
        df[col] = df[col].str.replace('\'', '')
        df[col] = df[col].str.replace(']', '')

    left_unique_in_cell(act_on_gene_oncogenic_gr_drug, 'Drugs', 'Drug')

    take_off_brackets_quotes(act_on_gene_oncogenic_gr_drug, 'Drug')

    # Group drug is a result. Merge it with act_res on Gene and Alteration, mode = left
    del act_res['Drugs']
    del act_res['Tumor Type']
    act_res.drop_duplicates(subset=['Gene', 'Alteration'], inplace=True)

    act_group_res = act_on_gene_oncogenic_gr_drug.merge(act_res, on=['Gene', 'Alteration'], left_index=True,
                                                        how='inner')  # 12 rows

    # *****FUSIONS ****

    allfusions = []
    with open('/mnt/volume_nyc3_01/output/GeneFuse_' + pat + '_result', 'r') as f:
        for i in f.readlines():
            if i.startswith("#Fusion:"):
                gs = i.split('Fusion: ')[1]
                genes2 = gs.split('___')
                ge1 = genes2[0].split('_ENST')[0]
                ge2 = genes2[1].split('_ENST')[0]
                fu = ge1 + '-' + ge2
                fu_reversed = ge2 + '-' + ge1
                # append splitted:
                allfusions.append(fu)
                allfusions.append(fu_reversed)

    # 4. Interpret fusions:

    act = pd.read_csv('db/oncokb_biomarker_drug_associations-2.tsv', sep='\t')

    # 4.1. Fusions with particular genes pair: contains - and Fusion

    act_fusion_concr = act[(act['Alterations'].str.contains('Fusion')) & (act['Alterations'].str.contains('-'))]

    search_list = '|'.join(allfusions)
    act_fusion_concr = act_fusion_concr[act_fusion_concr['Alterations'].str.contains(search_list)]

    # Add column with Alteration:

    act_fusion_concr = act_fusion_concr.assign(
        Alteration=act_fusion_concr['Alterations'].str.split(' ').str[0] + ' Fusion')

    # 4.2. Fusions with just all fusions participating this gene:
    act_all_fusions = act[act['Alterations'] == 'Fusions']

    act_all_fusions['Fus_flag'] = 0
    act_all_fusions['Alteration'] = ' '

    for pair in allfusions:
        ge1 = pair.split('-')[0]
        ge2 = pair.split('-')[1]
        act_all_fusions.loc[act_all_fusions['Gene'] == ge1, 'Fus_flag'] = 1
        act_all_fusions.loc[act_all_fusions['Gene'] == ge2, 'Fus_flag'] = 1
        # Add user alteration value:
        act_all_fusions.loc[act_all_fusions['Gene'] == ge1, 'Alteration'] = pair + ' Fusion'
        act_all_fusions.loc[act_all_fusions['Gene'] == ge2, 'Alteration'] = pair + ' Fusion'

    act_all_fusions = act_all_fusions.loc[act_all_fusions['Fus_flag'] == 1]
    del act_all_fusions['Fus_flag']

    # 4.3. Concat all 2 dfs about Fusions:

    fus = pd.concat([act_all_fusions, act_fusion_concr], sort=True)

    # Group:
    fus_gr_drug = fus.groupby(['Gene', 'Tumor Type', 'Alteration'])['Drugs'].apply(', '.join).reset_index()

    left_unique_in_cell(fus_gr_drug, 'Drugs', 'Drug')

    take_off_brackets_quotes(fus_gr_drug, 'Drug')

    # Group drug is a result. Merge it with act_res on Gene and Alteration, mode = left
    del fus['Drugs']
    del fus['Tumor Type']
    fus.drop_duplicates(subset=['Gene', 'Alteration'], inplace=True)

    fus_group_res = fus_gr_drug.merge(fus, on=['Gene', 'Alteration'], left_index=True, how='inner')  # 12 rows

    fus_group_res['FATHMM score'] = 1  # just for not empty. Must be int-float
    fus_group_res['FATHMM prediction'] = '-'
    fus_group_res['Pubmed_PMID'] = '-'

    # Concat with previous:

    act_group_res = pd.concat([act_group_res, fus_group_res], sort=True)

    # **** END Fusions insert ****

    # Cosmetics:
    '''
    'Gene', 'Tumor Type', 'Alteration', 'Drug', 'Alterations', 'Chromosome', 'FATHMM prediction', 'FATHMM score', 'Level', 
    'Mutation somatic status', 'Pubmed_PMID', 'Short_Annotation' '''

    # columns order
    act_group_res = act_group_res[
        ['Gene', 'Alteration', 'FATHMM prediction', 'FATHMM score', 'Alterations', 'Drug', 'Tumor Type', 'Level',
         'Pubmed_PMID', ]]

    # pubmed links
    act_group_res['Pubmed_PMID'] = act_group_res['Pubmed_PMID'].astype(str)
    act_group_res['Pubmed_PMID'] = act_group_res['Pubmed_PMID'].str.replace('\.0', '')

    act_group_res.reset_index(inplace=True)
    act_group_res.index = act_group_res.index + 1
    del act_group_res['index']

    act_group_res['FATHMM score'] = act_group_res['FATHMM score'].round(decimals=2)

    # Rename categories:
    new_evidences = {'1': 'A',
                     '2': 'B',
                     '2A': 'B',
                     '3A': 'C',
                     '3': 'C',
                     '4': 'D',
                     'R1': 'A-Resist',
                     'R2': 'C-Resist'}

    act_group_res['Level'].replace(new_evidences, inplace=True)

    act_group_res.rename(columns={'Alterations': 'Drug Indicators'}, inplace=True)
    act_group_res.rename(columns={'Pubmed_PMID': 'Pubmed PMID'}, inplace=True)

    # MSI
    # 1. MMR Genes ctitical mutations:

    msi_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2']

    cosmic_msi = cosmic_annotated_patho[cosmic_annotated_patho['Gene'].isin(msi_genes)]
    # or not
    # cosmic_msi = cosmic_annotated_patho[cosmic_annotated_patho['Gene'].isin(msi_genes)]

    if len(cosmic_msi.index) == 0:
        # make empty string with a message that didnt't found any:
        cosmic_msi_string = pd.DataFrame()

        cosmic_msi_string.at[len(act_group_res) + 1, 'Gene'] = 'MMR Genes (for MSI)'
        cosmic_msi_string['Alteration'] = 'Not found'
        cosmic_msi_string['FATHMM prediction'] = '-'
        cosmic_msi_string['FATHMM score'] = '-'
        cosmic_msi_string['Drug Indicators'] = 'Critical mutations'
        cosmic_msi_string['Drug'] = 'No immunotherapy'
        cosmic_msi_string['Tumor Type'] = 'All solid tumors'
        cosmic_msi_string['Level'] = 'A'
        cosmic_msi_string['Pubmed PMID'] = '28596308'

        overall = act_group_res.append(cosmic_msi_string)

    else:

        # Modify strings :

        cosmic_msi['Drug Indicators'] = 'MMR Genes Mutations'
        cosmic_msi['Drug'] = 'Pembrolizumab, Nivolumab, Nivolumab + Ipilimumab'
        cosmic_msi['Tumor Type'] = 'All solid tumors'
        cosmic_msi['Level'] = 'A'
        cosmic_msi['Pubmed PMID'] = '28596308'

        # New cols order - to get canonical:
        cosmic_msi = cosmic_msi[
            ['Gene', 'Alteration', 'FATHMM prediction', 'FATHMM score', 'Drug Indicators', 'Drug', 'Tumor Type',
             'Level', 'Pubmed PMID']]

        # Concat possible will be after all MSI and TMB filling:
        # overall = pd.concat([act_group_res, cosmic_msi], sort=True)
        overall = act_group_res.append(cosmic_msi)

        # Cosmetics:
        overall = overall[
            ['Gene', 'Alteration', 'FATHMM prediction', 'FATHMM score', 'Drug Indicators', 'Drug', 'Tumor Type',
             'Level', 'Pubmed PMID']]

        overall.reset_index(inplace=True)
        overall.index += 1
        del overall['index']

    # MSI:  2. Calculated from the script: insert to the result table

    p = pd.read_csv('/mnt/volume_nyc3_01/output/msisensor2_output_' + pat + '.tumor.prefix', sep='\t')

    '''Total_Number_of_Sites  Number_of_Somatic_Sites      %
    '''
    Total_Number_of_Sites = p.at[0, 'Total_Number_of_Sites']
    Number_of_Somatic_Sites = p.at[0, 'Number_of_Somatic_Sites']
    msi_percent = p.at[0, '%']

    # Empty string:
    msi_perc_string = pd.DataFrame()

    # Invariant strings:
    msi_perc_string.at[len(overall) + 1, 'Gene'] = 'MSI (Somatic sites)'
    msi_perc_string['FATHMM prediction'] = '-'
    msi_perc_string['FATHMM score'] = '-'
    msi_perc_string['Drug Indicators'] = '> 20.00 %'
    msi_perc_string['Tumor Type'] = 'All solid tumors'
    msi_perc_string['Level'] = 'A'
    msi_perc_string['Pubmed PMID'] = '28596308'

    # Condition:
    if msi_percent < 20:
        # make empty string with a message that lower:

        msi_perc_string['Alteration'] = 'Stable: ' + str(msi_percent) + '% (' + str(
            Number_of_Somatic_Sites) + '/' + str(Total_Number_of_Sites) + ')'
        msi_perc_string['Drug'] = 'No immunotherapy'

    else:

        msi_perc_string['Alteration'] = 'High: ' + str(msi_percent) + '% (' + str(Number_of_Somatic_Sites) + '/' + str(
            Total_Number_of_Sites) + ')'  # 2
        msi_perc_string['Drug'] = 'Pembrolizumab, Nivolumab, Nivolumab + Ipilimumab'  # 1

    overall2 = overall.append(msi_perc_string)

    # TMB

    # For TMB calc:
    tmb_nonsyn_total = len(cosmic_annotated.index)

    # App string to overall table:

    tmb_string = pd.DataFrame()

    tmb_string.at[len(overall2) + 1, 'Gene'] = 'TMB - Non-synonymous mutations'
    tmb_string['Alteration'] = str(tmb_nonsyn_total)
    tmb_string['FATHMM prediction'] = '-'
    tmb_string['FATHMM score'] = '-'
    tmb_string['Drug Indicators'] = '> 102 (KEYNOTE-012/028), > 100 (CheckMate 038)'
    tmb_string['Tumor Type'] = 'All solid tumors'
    tmb_string['Level'] = 'A'
    tmb_string['Pubmed PMID'] = '30664300'

    if tmb_nonsyn_total > 100:
        tmb_string['Drug'] = 'Pembrolizumab, Nivolumab + Ipilimumab'

    else:
        tmb_string['Drug'] = 'No immunotherapy'

    overall3 = overall2.append(tmb_string)

    # Finish - save result to folder:
    overall3.to_csv('/mnt/volume_nyc3_01/done/' + pat + '.csv')


def action_at_upload(value):
    """View uploaded patient on tab 2"""

    # Get calculated patient file:
    df = pd.read_csv('/mnt/volume_nyc3_01/done/' + str(value), index_col=0)
    msi_number = float(df.loc[df['Gene'] == 'MSI (Somatic sites)', 'Alteration'].iloc[0].split('Stable: ')[1].split('%')[0])
    tmb_number = int(df.loc[df['Gene'] == 'TMB - Non-synonymous mutations', 'Alteration'].iloc[0])

    # Condition - for different values:
    if 3 > 0:
        def output_for_done(value):
            return html.Div([
                html.Div([
                    html.H4(
                        className='what-is',
                        children='TUMOR COMPLEX GENOMIC PROFILING REPORT'
                    ),
                    html.H5(
                        """
                        PATIENT:
                        """
                    ),
                    html.P(
                        """                    
                        Name: Mikhail S. Yakovlev
                        """
                    ),
                    html.P(
                        """                    
                        Code: 2974-HJKD-20204389-
                        """ + str(value).split('.csv')[0]
                    ),
                    html.P(
                        """                    
                        Disease: Colorectal cancer
                        """
                    ),
                    html.P(
                        """                    
                        Date of Birth: 11.06.63      
                        """
                    ),
                    html.H5(
                        """
                        MATERIALS AND METHODS:
                        """
                    ),
                    html.P(
                        """
                        The biosample (FFPE paraffin block) was sequenced using NGS paired-end method on Illumina NovaSeq 7000, the enrichment was performed using the Agilent v6 kit. Alignment to the GRCh38 / hg38 reference genome was performed using Mutect2 from GATK software package by Broad Institute.
                        Amino acid changes were detected using Variant Effect Predictor of ENSEMBL Project by EBBL-EBI. Pathogenic variants were detected using COSMIC by Sanger Institute, UK and FATHMM by University of Bristol, UK databases [4]. Actionable variants were detected using clinicaltrials.gov database.
                        MSI was calculated using MSIsensor2 software as recommended by Standing Operating Procedure MCCRD-SOP0057 (EXT) by National Cancer Institute [2,7].
                        Fusions were detected using GeneFuse software [7].
                        """
                    ),
                    html.H5(
                        """
                        ACTIONABLE VARIANTS DETECTED:
                        """
                    ),
                ]),
                html.Div([
                    dash_table.DataTable(
                        id='datatable-interactivity',
                        columns=[
                            {"name": i, "id": i, "deletable": True} for i in df.columns
                        ],
                        style_cell={'textAlign': 'left',
                                    'whiteSpace': 'normal',
                                    'height': 'auto',
                                    },
                        data=df.to_dict('records'),
                        editable=True,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        column_selectable="single",
                        # row_selectable="multi",
                        row_deletable=True,
                        selected_columns=[],
                        selected_rows=[],
                        page_action="native",
                        page_current=0,
                        # page_size=10,

                    )

                ]),

                html.Div([
                    html.H5(
                        children='CONCLUSION:'
                    ),
                    html.P(
                        """
                    Detected alterations can be neutralized with the appropriate drugs indicated in the table in the "off-label" mode (the diagnosis for which the application was tested does not correspond to the patient's diagnosis)

                        """
                    ),



                    html.H5(
                        children='IMMUNE THERAPY:'
                    ),

                    html.P(
                        """
                    Detected alterations can be neutralized with the appropriate drugs indicated in the table in the "off-label" mode (the diagnosis for which the application was tested does not correspond to the patient's diagnosis)
        
                        """ + ' efe'  + str(msi_number) + ' elfgk' + str(tmb_number)

                        # Get MSI percent:

                    ),




                    html.H5(
                        """
                        REFERENCES:
                        """
                    ),
                    html.P(
                        """
                            1.        Rizvi,N.A.et al. Cancer immunology. Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer. Science, 348(6230):124- 128, April 2015."""),
                    html.P(
                        """
                2.        Salipante,S.J.et al. Microsatellite Instability Detection by Next Generation Sequencing. Clinical Chemistry, February 2014. """),
                    html.P(
                        """
                3.        D. T. Le, et al. Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. Science. Published Online 8 June 2017. DOI: 10.1126/science.aan6733. """),
                    html.P(
                        """
                4.        COSMIC: somatic cancer genetics at high-resolution. Simon A. Forbes David Beare Harry Boutselakis Sally Bamford Nidhi BindalJohn Tate Charlotte G. Cole Sari Ward Elisabeth Dawson Laura Ponting. Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D777–D783,https://doi.org/10.1093/nar/gkw1121"""),
                    html.P(
                        """
                5.        RNA sequence analysis reveals macroscopic somatic clonal expansion across normal tissues Yizhak K, Aguet F, Kim J, Hess JM, Kübler K et al. Science. 07 June 2019. 364(6444). doi:10.1126/science.aaw0726"""),
                    html.P(
                        """
                6.        A Pan-Cancer Analysis of Enhancer Expression in Nearly 9000 Patient Samples, Han Chen, Chunyan Li 4, Xinxin Peng, John N. Weinstein. Cell, VOLUME 173, ISSUE 2, P386-399.E12, APRIL 05, 2018"""),
                    html.P(
                        """
                7.        GeneFuse: detection and visualization of target gene fusions from DNA sequencing data. Shifu Chen, Ming Liu, Tanxiao Huang, Wenting Liao, Mingyan Xu, Jia Gu. Int J Biol Sci. 2018; 14(8): 843–848. Published online 2018 May 22. doi: 10.7150/ijbs.24626
                8.        https://pdmr.cancer.gov/content/docs/MCCRD_SOP0057_MSI_Status_from_WES.pdf
            """
                    ),
                    html.H5(
                        """
                        CONTACTS:
                        """
                    ),
                    html.P(
                        """
                            Electronically signed by Dmitrii K. Chebanov, CEO

                            Read more about us here:
                            http://bioalg.com/
                            Contact us at:
                            info@bioalg.com
                        """
                    ),
                ]),

            ])

    return output_for_done(value)


class Semaphore:
    def __init__(self, filename='semaphore.txt'):
        self.filename = filename
        with open(self.filename, 'w') as f:
            f.write('done')

    # With file name which in progress:
    # def lock(self, value):
    # with open(self.filename, 'w') as f:
    # f.write('working ' + str(value))

    def lock(self):
        with open(self.filename, 'w') as f:
            f.write('working')

    def unlock(self):
        with open(self.filename, 'w') as f:
            f.write('done')

    def is_locked(self):
        return open(self.filename, 'r').read() == 'working'


semaphore = Semaphore()


# 1. Self-check
@app.callback(
    Output('status-message', 'children'),
    [Input('interval', 'n_intervals'),
     ],
    # prevent_initial_call=True
)
def self_check(n):
    return 'Running...' if semaphore.is_locked() else 'Provide the data and press "Calculate"'


# 2. Button press for message (maybe not need, delete ):
@app.callback(
    Output('err-message', 'children'),
    [Input('input-seq-1', 'value'),
     Input('input-seq-2', 'value'),
     Input('calculate', 'n_clicks'),
     # This status message turns off (wakes) staying non-actual message
     Input('status-message', 'children')
     ],
    [State('input-name', 'value')],
    prevent_initial_call=True
)
def status_show(lin1, lin2, n_clicks, m, value):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'calculate' in changed_id:
        if (lin1 == None) | (lin2 == None) | (value == None):
            status = 'Please fill all the fields!'
            return status
        else:
            if semaphore.is_locked():
                status = 'Resource is busy. Try later'
                return status
            else:
                # action_at_click(value)
                status = 'Successfully submitted ' + str(value) + '. Find it on the "Calculated" tab in 6-8 hours.'
                return status


# 3. Button press for action (maybe not need, delete - combine with previous):
@app.callback(
    Output('result-print', 'children'),
    [Input('calculate', 'n_clicks'),
     Input('input-seq-1', 'value'),
     Input('input-seq-2', 'value')
     ],
    [State('input-name', 'value')]
)
def calculation(n_clicks, seq1, seq2, value):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'calculate' in changed_id:
        if (seq1 == None) | (seq2 == None) | (value == None):
            raise PreventUpdate
        else:
            if not semaphore.is_locked():
                action_at_click(seq1, seq2, value)


# 4. Scroll the calculated patients:
@app.callback(
    Output('result-field', 'children'),
    [Input('calculated-dropdown', 'value')
     ],
    prevent_initial_call=True   #or not???
)
def calculated_file_show(patient):
    #action_at_upload(patient)
    return action_at_upload(patient)




@app.callback(
    Output('callback-output', 'children'),
    [Input('alignment-file-upload', 'isCompleted')],
    [State('alignment-file-upload', 'fileNames'),
     State('alignment-file-upload', 'upload_id')],
)
def callback_on_completion(iscompleted, filenames, upload_id):
    if not iscompleted:
        return

    out = []
    if filenames is not None:
        if upload_id:
            root_folder = Path(UPLOAD_FOLDER_ROOT) / upload_id
        else:
            root_folder = Path(UPLOAD_FOLDER_ROOT)

        for filename in filenames:
            file = root_folder / filename
            out.append(file)
        return html.Ul([html.Li(str(x)) for x in out])

    return html.Div("No Files Uploaded Yet!")

# Handle event data
@app.callback(
    Output("alignment-events", "children"),
    [Input("alignment-chart", "eventDatum")]
)
def event_data_select(data):
    if data is None:
        data = '{}'

    data = json.loads(data)

    if len(data.keys()) == 0:
        return 'No event data to display.'

    return [
        html.Div('- {}: {}'.format(key, data[key]))
        for key in data.keys()
    ]


if __name__ == '__main__':
    app.run_server(
        port=80,  # 80 - default for web-browser, 8050 - default for dash
        host='138.197.1.131',
        debug=True
    )