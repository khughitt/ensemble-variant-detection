"""
Read-mapping classes
"""
class Mapper(object):
    """Base read mapper class"""
    def __init__(self):
        """Create a new Mapper instance"""
        pass

class BWAMapper(Mapper):
    """Burrows-Wheeler Aligner Mapper class"""
    def __init__(self):
        pass

    def run(self):
        pass
