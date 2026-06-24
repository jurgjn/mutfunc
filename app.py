
"""
uv run app.py
uv run gunicorn app:server --bind 0.0.0.0:8050
"""
from pprint import pprint
import os
import pandas as pd
import dash
from dash import Dash, html, Input, Output
import dash_ag_grid as dag
import datetime
import dash_molstar
from dash_molstar.utils import molstar_helper

import dash_bootstrap_components as dbc

from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

#from mutationimpact_template import register_template, PREDICTION_COLORS

#register_template(set_default=True)

MUTFUNC_BUILD = os.environ.get("MUTFUNC_BUILD", f'{datetime.datetime.now().strftime("%y.%m.%d")}-dev')

footer_ = html.Footer(
    dbc.Container(
        dbc.Row([
            dbc.Col(
                dcc.Markdown(
                    "Prof. Dr. Pedro Beltrao  \n"
                    "[E-Mail](mailto:beltrao@imsb.biol.ethz.ch)  \n"
                    "Otto-Stern-Weg 3, 8093 Zurich, Switzerland",
                    link_target="_blank",
                ),
            ),
            dbc.Col(
                html.Span(MUTFUNC_BUILD),
                className="text-end",
            ),
        ]),
        fluid=True,
        style={"padding": "1rem"},
    ),
    style={
        "marginTop": "2rem",
        "borderTop": "1px solid #dee2e6",
        "fontSize": "0.875rem",
        "color": "#6c757d",
    },
)

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], use_pages=True)
app.title = "mutfunc"
app.layout = html.Div([
    dcc.Store(id="variant-list", storage_type="session"),
    dbc.Navbar(
        dbc.Container([
            dbc.NavbarBrand("mutfunc - precomputed mechanistic consequences of missense mutations", href="/"),
            dbc.Nav([
                dbc.NavLink("Home", href="/"),
                dbc.NavLink("About", href="/about"),
                ],
                navbar=True,
                className="ms-auto",  # pushes nav links to the right
            ),
            ],
            fluid=True,  # change to fluid=True so it spans full width
        ),
        className="bg-canvas",
        sticky="top",
        dark=True,
    ),
    html.Div(
        dash.page_container,
        style={
            "maxWidth": "1500px",
            "margin": "0 auto",
            "padding": "0 1rem",
        }
    ),
    footer_,
])

server = app.server

if __name__ == '__main__':
    app.run(debug=True)
