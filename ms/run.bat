ECHO
SET mydoc=alaska_ms_master
pdflatex "%mydoc%.tex"
pandoc --mathjax --reference-docx=c:\~\localtex\reference.docx "%mydoc%.tex" -o "alaska_ms_kfj.docx"

    for %%f in (*.aux) do (
            DEL "%%f"
    )
    for %%f in (*.log) do (
            DEL "%%f"
    )
    for %%f in (*_latexmk) do (
            DEL "%%f"
    )
    for %%f in (*.fls) do (
            DEL "%%f"
    )
alaska_ms_kfj.docx

exit