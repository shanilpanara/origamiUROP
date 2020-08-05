from sandbox.readers import ReaderTopConf
from origamiUROP.oxdna import System
import filecmp

def test_oxDNA_reader():
    conf = "./oxdna.out.conf"
    top = "./oxdna.out.top"
    system = ReaderTopConf(top, conf).return_system()
    system.write_oxDNA("reader")
    # this doesn't quite work just yet
    assert filecmp.cmp("./oxdna.out.conf", "./oxdna.reader.conf")

if __name__ == "__main__":
    def test_oxDNA_reader()