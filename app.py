
"""
uv run app.py
uv run gunicorn app:server --bind 0.0.0.0:8050
"""
from pprint import pprint

import pandas as pd
import dash
from dash import Dash, html, Input, Output
import dash_ag_grid as dag

import dash_molstar
from dash_molstar.utils import molstar_helper

import dash_bootstrap_components as dbc

from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from mutationimpact_template import register_template, PREDICTION_COLORS

register_template(set_default=True)

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], use_pages=True)
app.title = "mutfunc"
app.layout = html.Div([
    dcc.Store(id="variant-list", storage_type="session"),
    dash.page_container,
])

# Force Bootstrap's dark color mode at the <html> level.
app.index_string = """<!DOCTYPE html>
<html lang="en" data-bs-theme="dark">
<head>{%metas%}<title>{%title%}</title>{%favicon%}{%css%}</head>
<body>
  {%app_entry%}
  <footer>{%config%}{%scripts%}{%renderer%}</footer>
</body>
</html>"""

'''
app = Dash(external_stylesheets=[
    #dbc.themes.SUPERHERO
    "/assets/mutfunc.css",
], use_pages=True)
'''

server = app.server

if __name__ == '__main__':
    app.run(debug=True)
