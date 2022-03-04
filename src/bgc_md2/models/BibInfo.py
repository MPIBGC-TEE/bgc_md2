import inspect



class BibInfo():
    def __init__(self,**kwargs):
        fields = {
            'name',
            'longName',
            'version',
            'entryAuthor',
            'entryAuthorOrcid',
            'entryCreationDate',
            'doi',
            'further_references',
            'sym_dict',
            'func_dict',
            'exprs_dict'
        }
        keyset = set(kwargs.keys())
        if not(keyset.issubset(fields)):
            raise Exception("The fields "+str(keyset.difference(fields))+" are not allowed")

        for k, v in kwargs.items():
            setattr(self, k, v)
