from origamiUROP.readers import ReaderTopConf
from origamiUROP.oxdna import System
import filecmp
from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def test_oxDNA_reader():
    top_input = f"{ROOT}/oxdna.test.top"
    conf_input = f"{ROOT}/oxdna.test.conf"
    print(top_input, conf_input)
    system = ReaderTopConf(top_input, conf_input).return_system()

    system.write_oxDNA("reader", ROOT)

    top_out = f"{ROOT}/oxdna.reader.top"
    conf_out = f"{ROOT}/oxdna.reader.conf"
    assert filecmp.cmp(top_input, top_out), "topology files don't match"
    assert filecmp.cmp(conf_input, conf_out), "conf files don't match"


if __name__ == "__main__":
    test_oxDNA_reader()
