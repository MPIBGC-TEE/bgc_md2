from IPython.display import HTML

display(HTML("<style>.container { width:100% !important; }</style>"))

# %load_ext autoreload
# %autoreload 2
import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()

model_list = h.ModelListGridBox(
    inspection_box=model_inspection,
    explicit_exclude_models=frozenset({'CARDAMOM'})
)
model_list

model_inspection

model_list.inspect_mvs(model_list.names[0])
