from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import fitter
import plotly.graph_objects as go


FOOTER_STYLE = {
    "position": "fixed",
    "bottom": 0,
    "left": 0,
    "right": 0,
    "height": 50,
    "padding": "1rem 1rem",
#    "background-color": "white",
}

class PSDAnalyzerUI:
    def __init__(self):
        self.layout = html.Div([
            dbc.Row([
                dbc.Col([html.H1("Particule Size Distribution Analyzer", style={'textAlign':'center'})]),
                dbc.Col([
                    dbc.Button("Get help!", href="https://docs.hydrogeotools.com/psd-analyzer.html"),
                ], width="auto"),
            ], justify="end", align="center"),
            html.Hr(),
            dbc.Row([
                dbc.Col(self._visualization_pannel(), width=7),
                dbc.Col([
                    dbc.Row(self._calibration_pannel()),
                    dbc.Row(self._results_pannel())
                ]),
            ]),
            dbc.Row([self._footer()]),
        ])
        return
    
    def _visualization_pannel(self):
        header = [
            html.H3("Experimental data", className="card-title"),
        ]
        upload = [
            html.H6("Upload your particule size distribution: ", id="upload-text"),
            dcc.Upload(
                id='upload-data',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select a file'),
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
            ),
        ]
        graph = [
            html.H6("Visualize your data:", id="data-vis-text"),
            dcc.Graph(id='graph-content'),
            dbc.Button('Download a sample curve', id='button-sample-wrc', n_clicks=0),
            dcc.Download(id="download-sample-curve"),
        ]
        visualization_card = dbc.Card(
            dbc.CardBody(
                header +
                upload +
                graph
            ),
            style={"width": "100%"}
        )
        return visualization_card
    
    def _calibration_pannel(self):
        header = [
            html.H3("Calibration properties",className="card-title")
        ]
        fit_type = [
            html.H6("Select fit type: ", id="fit-text"),
            dbc.RadioItems(["Best fit with RMSE (deterministic)", "Quantile regression (statistical)"], "Best fit with RMSE (deterministic)", id="radio-fit-selector"),
            html.Hr(),
        ]
        opt_button = [dbc.Button('Optimize', id='optimize_button', n_clicks=0)]
        calibration_card = dbc.Card(
            dbc.CardBody(
                header +
                fit_type +
                opt_button
            ),
            style={"width": "100%"}
        )
        return calibration_card
    
    def _results_pannel(self):
        header = [
            html.H3("Results"),
        ]
        tab_fitted_params = [
            #html.H6("Calibrated Fredlund PSD parameters:", id="out-fitted-params-text"),
            html.Div(["Please optimize to view the fitted parameters."], id='fitted-result-data'),
        ]
        tab_geo_params = [
            #html.H6("Derived geotechnical parameters:", id="out-geo-params-text"),
            html.Div(["Please optimize to view the derived geotechnical parameters."], id='geo-result-data'),
        ]
        tab_hydro_params = [
            html.H6("Saturated hydraulic conductivity", id="out-hydro-params-text"),
            html.P("Prediction of the saturated hydraulic conductivity following Mbonimpa et al. (2002) considering a granular soil (i.e. non-plastic and non-cohesive). To recover hydraulic conductivity, put permeability in m^2 and multiply by 10^7."),
            dcc.Graph(id='graph-hydro-content')
        ]
        tabs = [dbc.Tabs([
            dbc.Tab(tab_fitted_params, label="Calibrated PSD parameters"),
            dbc.Tab(tab_geo_params, label="Derived geotechnical data"),
            dbc.Tab(tab_hydro_params, label="Derived hydrogeological data"),
        ])]
        download = [
            html.H6("Download curve(s) in output format:", id="out-results-format-text"),
            dcc.Dropdown([
                "Excel (.xlsx)",
                "CSV (.csv)",
            ], "CSV (.csv)", id="dropdown-out-format-selector"),
            dbc.Button('Download', id='download-results-btn', n_clicks=0),
            dcc.Download(id="download-results"),
        ]
        results_card = dbc.Card(
            dbc.CardBody(
                header
                + tabs
                + download
            )
        )
        return results_card
    
    def _footer(self):
        contents = [ html.P([
            "Application created by HydroGeoTools. See ", 
            html.A("www.hydrogeotools.com", href="https://www.hydrogeotools.com"),
            " for more contents.",
        ]) ]
        footer = html.Footer(contents, style=FOOTER_STYLE)
        return footer
    
    def packLayout(self):
        return dbc.Container(self.layout, fluid=True)
        




def _present_psd_params(fit_type, res):
    if fit_type == "Best fit with RMSE (deterministic)":
        fitted_table_header = html.Thead(
            html.Tr([
                html.Th("Fredlund PSD Parameters"),
                html.Th("Best fit (RMSE)")
            ])
        )
        fitted_table_body = html.Tbody([
            html.Tr([html.Td(name), html.Td(f"{10**res[0].x[i]:.6e}")]) for i,name in enumerate(fitter.model_output_name)
        ])
    elif fit_type == "Quantile regression (statistical)":
        fitted_table_header = html.Thead(
            html.Tr([
                html.Th("Fredlund PSD Parameters"),
                html.Th("5th Quantile"),
                html.Th("50th Quantile"),
                html.Th("95th Quantile")
            ])
        )
        resmin, res50, resmax = res
        fitted_table_body = html.Tbody([
            html.Tr(
                [html.Td(name), html.Td(f"{10**resmin.x[i]:.6e}"), html.Td(f"{10**res50.x[i]:.6e}"), html.Td(f"{10**resmax.x[i]:.6e}")]) for i,name in enumerate(fitter.model_output_name
            )
        ])
    children_fitted_params = html.Div([
        dbc.Table([fitted_table_header, fitted_table_body], bordered=True), 
        html.P("From the finner (5th quantile) to the coarser soil (95th quantile). The above parameters table is copy-pastable in Excel")
    ])
    return children_fitted_params


def _present_geo_params(fit_type, geo_results):
    if fit_type == "Best fit with RMSE (deterministic)":
        geo_table_header = html.Thead(
            html.Tr([
                html.Th("Geotechnical Parameters"),
                html.Th("Best fit (RMSE)")
            ])
        )
    elif fit_type == "Quantile regression (statistical)":
        geo_table_header = html.Thead(
            html.Tr([
                html.Th("Geotechnical Parameters"),
                html.Th("5th Quantile"),
                html.Th("50th Quantile"),
                html.Th("95th Quantile")
            ])
        )
    geo_table_body = html.Tbody([
        html.Tr([html.Td(name)] + [html.Td(f"{val:.2e}") for val in array]) for name,array in geo_results.items()
    ])
    children_geo_params = html.Div([
        dbc.Table([geo_table_header, geo_table_body], bordered=True), 
        html.P("From the finner (5th quantile) to the coarser soil (95th quantile). The above parameters table is copy-pastable in Excel")
    ])
    return children_geo_params
