import qrcode

def make_qr(text,filename):
    qrcode.make(text).save(filename)

for tup in [
        (
            'https://github.com/MPIBGC-TEE/bgc_md2#readme',
            "qr_github.pdf"
        ),
        (
            'https://github.com/MPIBGC-TEE/bgc_md2/blob/test/notebooks/illustrativeExamples/createModel.ipynb',
            'qr_tutorial.pdf'
        ),    
        (
            'mailto:markus.mueller.1.g@googlemail.com',
            'qr_mm.pdf'
        )    
    ]:
    make_qr(*tup)
