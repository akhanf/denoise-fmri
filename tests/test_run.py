import subprocess

def test_bids_01():
    subprocess.check_call(['denoise_fmri','tests/bids_01','tests/out_01','participant','-np'])
