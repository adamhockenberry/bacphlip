import pytest
import bacphlip

import os
import filecmp

def test_example_files():
    import pkg_resources
    file_dict = {}
    for example_file in ['genome_example.fasta',\
                         'genome_example.fasta.6frame',\
                         'genome_example.fasta.hmmsearch',\
                         'genome_example.fasta.hmmsearch.tsv',\
                         'genome_example.fasta.bacphlip']:
        example_path = pkg_resources.resource_filename('bacphlip', 'data/example_data/{}'.format(example_file))
        assert os.path.exists(example_path)    
        assert os.path.getsize(example_path) != 0
        file_dict[example_file] = example_path
    return file_dict

def test_support_files():
    assert os.path.exists(bacphlip.HMMER_DB)    
    assert os.path.exists(bacphlip.SKLEARN_CLASSIFIER)    
    assert os.path.getsize(bacphlip.HMMER_DB) != 0    
    assert os.path.getsize(bacphlip.SKLEARN_CLASSIFIER) != 0

def check_non_existent_file():
    import pkg_resources
    invalid_file = 'genome_example.fastaaaaaaa'
    example_path = pkg_resources.resource_filename('bacphlip', 'data/example_data/{}'.format(invalid_file))
    with pytest.raises(Exception):
        bacphlip.check_existing_file(example_path)

def test_translate(tmp_path):
    file_dict = test_example_files()
    tmp_six_frame_path = tmp_path / 'six_frame.fasta' 
    assert os.path.exists(tmp_six_frame_path) == False
    bacphlip.six_frame_translate(file_dict['genome_example.fasta'], tmp_six_frame_path)
    assert (os.path.exists(tmp_six_frame_path) == True) and (os.path.getsize(tmp_six_frame_path) != 0)
    assert filecmp.cmp(file_dict['genome_example.fasta.6frame'], tmp_six_frame_path, shallow=False)

def test_hmmsearch(tmp_path):
    file_dict = test_example_files()
    tmp_hmmsearch_out = tmp_path / 'tempy.fasta.hmmsearch' 
    assert os.path.exists(tmp_hmmsearch_out) == False
    bacphlip.hmmsearch_py(file_dict['genome_example.fasta.6frame'], tmp_hmmsearch_out)
    assert (os.path.exists(tmp_hmmsearch_out) == True) and (os.path.getsize(tmp_hmmsearch_out) != 0)
    with open(file_dict['genome_example.fasta.hmmsearch'], 'r') as temp:
        lines_a = len(temp.readlines()) 
    with open(tmp_hmmsearch_out, 'r') as temp:
        lines_b = len(temp.readlines()) 
    assert lines_a == lines_b

def test_process_hmmsearch(tmp_path):
    file_dict = test_example_files()
    tmp_hmmsearchtsv_path = tmp_path / 'six_frame.fasta.hmmsearch.tsv' 
    assert os.path.exists(tmp_hmmsearchtsv_path) == False
    bacphlip.process_hmmsearch(file_dict['genome_example.fasta.hmmsearch'], tmp_hmmsearchtsv_path)
    assert (os.path.exists(tmp_hmmsearchtsv_path) == True) and (os.path.getsize(tmp_hmmsearchtsv_path) != 0)
    assert filecmp.cmp(file_dict['genome_example.fasta.hmmsearch.tsv'], tmp_hmmsearchtsv_path, shallow=False)

def test_bacphlip(tmp_path):
    file_dict = test_example_files()
    tmp_bacphlip_path = tmp_path / 'six_frame.fasta.bacphlip' 
    assert os.path.exists(tmp_bacphlip_path) == False
    bacphlip.predict_lifestyle(file_dict['genome_example.fasta.hmmsearch.tsv'], tmp_bacphlip_path)
    assert (os.path.exists(tmp_bacphlip_path) == True) and (os.path.getsize(tmp_bacphlip_path) != 0)
    assert filecmp.cmp(file_dict['genome_example.fasta.bacphlip'], tmp_bacphlip_path, shallow=False)

def test_no_overwrite(tmp_path):
    file_dict = test_example_files()
    tmp_six_frame_path = tmp_path / 'six_frame.fasta' 
    assert os.path.exists(tmp_six_frame_path) == False
    bacphlip.six_frame_translate(file_dict['genome_example.fasta'], tmp_six_frame_path)
    assert (os.path.exists(tmp_six_frame_path) == True) and (os.path.getsize(tmp_six_frame_path) != 0)
    with pytest.raises(Exception):
        bacphlip.six_frame_translate(file_dict['genome_example.fasta'], tmp_six_frame_path)
