import argparse

def get_parser():
     
    parser = argparse.ArgumentParser(description='Hippocampal AutoTop BIDS App')
    parser.add_argument('bids_dir', help='The directory with the input dataset '
                        'formatted according to the BIDS standard.')
    parser.add_argument('output_dir', help='The directory where the output files '
                        'should be stored. If you are running group level analysis '
                        'this folder should be prepopulated with the results of the'
                        'participant level analysis.')
    parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                        'participant: Hippocampal unfolding and segmentation pipeline, '
                        'group: Reports and statistics (not implemented yet)',
                        choices=['participant', 'group'])
    parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                       'corresponds to sub-<participant_label> from the BIDS spec '
                       '(so it does not include "sub-"). If this parameter is not '
                       'provided all subjects should be analyzed. Multiple '
                       'participants can be specified with a space separated list.',
                       nargs="+")
    parser.add_argument('--session', help='Use the specified session, otherwise will process all sessions')

    #maybe remove suffix -- is not really typical to have a choice of suffix in an app (i.e. designed for particular suffixes)
#    parser.add_argument('--suffix', help='Only use images with the specified suffix entity in the filename',
#                        default=['T2w'],nargs='+')
    parser.add_argument('--acq', help='Only use images with the specified acq entity in the filename')
    parser.add_argument('--run', help='Only use images with the specified run entity in the filename')
    parser.add_argument('--search', help='Wildcard search term to locate in image filename. Use this option '
                        'when multiple images match for a subject')

    parser.add_argument('--skip_bids_validator', help='Whether or not to perform BIDS dataset validation',
                       action='store_true', default=True)


    return parser


