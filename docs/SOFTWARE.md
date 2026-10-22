# Software and local environment

R execution logs record R 4.5.2 on Windows. The versions below were detected in the local R library during source packaging; they are not a dependency lock for every historical run. Original session logs remain in the private project.

| R package | Detected version |
| --- | --- |
| base | 4.5.2 |
| data.table | 1.18.0 |
| digest | 0.6.39 |
| gbmt | 0.1.4 |
| ggplot2 | 4.0.1 |
| grDevices | 4.5.2 |
| graphics | 4.5.2 |
| grid | 4.5.2 |
| jsonlite | 2.0.0 |
| mice | 3.19.0 |
| mitml | 0.4-5 |
| parallel | 4.5.2 |
| patchwork | 1.3.2 |
| ragg | 1.5.0 |
| sandwich | 3.1-1 |
| stats | 4.5.2 |
| survival | 3.8-3 |
| tools | 4.5.2 |
| utils | 4.5.2 |
| zip | 2.3.3 |

Python imports are listed in `requirements.txt`. The following versions were detected in the packaging runtime; some earlier presentation scripts used a different local or vendored runtime.

| Python distribution | Detected version |
| --- | --- |
| numpy | 2.3.5 |
| pandas | 3.0.1 |
| Pillow | 12.3.0 |
| python-docx | 1.2.0 |
| lxml | 6.1.1 |
| reportlab | 4.4.9 |
| pypdfium2 | 5.14.0 |
| resvg-py | Not installed in the packaging runtime; required only by relevant presentation helper |

Some original renderers use Windows Arial/Calibri font paths; adapt these on other systems. PowerShell helpers support Windows CPU/process handling. Word-dependent formatting helpers are presentation utilities, not statistical dependencies. PostgreSQL schemas and versions must match the locally licensed database installation.

Clinical extraction, model fitting and data-dependent analyses were not executed for packaging. Python/SQL/R static parsing does not replace running the pipeline with authorised local data.
