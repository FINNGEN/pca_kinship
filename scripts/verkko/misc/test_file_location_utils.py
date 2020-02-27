import file_location_utils as flu
import os

def test_flu():
    p = flu.get_dir_to_cur_module(__file__)
    assert os.path.exists(p)
    assert flu.get_dir_to_cur_module(__file__)[0] == "/"
 

if __name__ == "__main__":
    test_flu()
