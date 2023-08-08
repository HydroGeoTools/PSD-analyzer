import base64
import datetime
import io

from dash import Dash, html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import pandas as pd
import numpy as np

import fitter
import main_ui as UI_file



app = Dash(__name__, title="PSD-Analyzer", external_stylesheets=[dbc.themes.BOOTSTRAP])
main_ui = UI_file.PSDAnalyzerUI()
app.layout = main_ui.packLayout()


def parse_contents(contents, filename):
    if contents is None: return
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')),
                sep=None,
                engine="python",
            )
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
        #clean data
        mask_nan = df.isna().sum(axis=1) == 0
        df = df[mask_nan]
    except Exception as e:
        print(e)

    return df


@app.callback(
    Output('graph-content', 'figure', allow_duplicate=True),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    prevent_initial_call=True
)
def update_graph(contents,filename):
    data = parse_contents(contents, filename)
    header = data.columns.tolist()
    xdata = np.array(data.iloc[:,0])
    ydata = np.array(data.iloc[:,1])
    fig = go.Figure(data=go.Scatter(x=xdata, y=ydata, mode='markers'))
    fig.update_xaxes(type="log")
    fig.update_layout(xaxis_title=header[0], yaxis_title=header[1], legend={"xanchor":"left", "yanchor":"top"})
    return fig


@app.callback(
    Output('graph-content', 'figure'),
    Output('fitted-result-data', 'children'),
    Output("geo-result-data", "children"),
    Output("graph-hydro-content", "figure"),
    Input('optimize_button', 'n_clicks'),
    State('upload-data', 'contents'),
    State('upload-data', 'filename'),
    State('radio-fit-selector', 'value'),
    #State('btn-download-res', 'hidden'),
    prevent_initial_call=True,
)
def optimize(btn, contents, filename, fit_type):
    if contents is None: 
        return go.Figure(), html.Div()
    
    #parse inputs
    data = parse_contents(contents, filename)
    header = data.columns.tolist()
    xdata = np.array(data.iloc[:,0])
    ydata = np.array(data.iloc[:,1])
    
    #fit
    if fit_type == "Best fit with RMSE (deterministic)":
        res = [fitter.fit(xdata, ydata)]
        loss = np.sqrt(res[0].fun)
        x_th = np.logspace(res[0].x[4], np.log10(np.max(xdata)), 200)
        x_th = np.append(x_th, xdata)
        x_th.sort()
        graph_data = [
            {
                "x" : x_th,
                "y" : fitter.Fredlund_PSD(x_th, *res[0].x),
                "line" : {"dash":"solid", "color":"red"}, 
                "marker" : None,
                "name" : "Fitted Fredlund PSD (Best fit - RMSE)",
            }
        ]
        
    elif fit_type == "Quantile regression (statistical)":
        res50 = fitter.fit_quantile(xdata, ydata, 0.50)
        resmax = fitter.fit_quantile(xdata, ydata, 0.95, res50.x)
        resmin = fitter.fit_quantile(xdata, ydata, 0.05, res50.x)
        res = [resmin, res50, resmax]
        x_th = np.logspace(resmax.x[4], np.log10(np.max(xdata)), 200)
        x_th = np.append(x_th, xdata)
        x_th.sort()
        loss = np.sqrt(resmax.fun**2 + res50.fun**2 + resmin.fun**2)
        graph_data = [
            {
                "x" : x_th,
                "y" : fitter.Fredlund_PSD(x_th, *resmin.x),
                "line" : {"dash":"solid", "color":"silver"}, 
                "marker" : None,
                "name" : f"5th Quantile",
            },
            {
                "x" : x_th,
                "y" : fitter.Fredlund_PSD(x_th, *res50.x),
                "line" : {"dash":"solid", "color":"red"}, 
                "marker" : None,
                "name" : f"50th Quantile",
            },
            {
                "x" : x_th,
                "y" : fitter.Fredlund_PSD(x_th, *resmax.x),
                "line" : {"dash":"solid", "color":"gold"}, 
                "marker" : None,
                "name" : f"95th Quantile",
            }
        ]
    
    # print text results
    children_fitted_params = UI_file._present_psd_params(fit_type, res)
    
    #compute D_10, D_60, Cu and Cc
    D_10 = np.array([fitter.Dx(10,fitter.Fredlund_PSD, r.x) for r in res])
    D_30 = np.array([fitter.Dx(30,fitter.Fredlund_PSD, r.x) for r in res])
    D_60 = np.array([fitter.Dx(60,fitter.Fredlund_PSD, r.x) for r in res])
    geo_results = {
        "D_10" : D_10,
        "D_30" : D_30,
        "D_60" : D_60,
        "Cu" : D_60/D_10,
        "Cc" : D_30*D_30/(D_60*D_10),
    }
    children_geo_params = UI_file._present_geo_params(fit_type, geo_results)
    
    #predict saturated perm following eq. 13 in Mbonimpa et al. (2002)
    #D_10 must be in cm!
    porosity = np.linspace(0.1, 0.6, 100)
    if fit_type == "Best fit with RMSE (deterministic)":
        labels = ["Best fit"]
    elif fit_type == "Quantile regression (statistical)":
        labels = ["Coarser soil", "Median soil", "Finner soil"]
    K_g = {
        label : fitter.predict_sat_perm_granular(porosity, geo_results["Cu"][i], D_10[i]) for i,label in enumerate(labels)
    }
    color = {"Finner soil":"gold","Median soil":"red","Coarser soil":"silver", "Best fit":"red"}
    
    fig_hydro = go.Figure(
        data=[
            go.Scatter(
                x=porosity,
                y=the_K_g,
                line={"dash":"solid", "color":color[label]},
                marker=None,
                name=label
            ) for label,the_K_g in K_g.items()
        ],
        layout = {
            "legend": {"xanchor":"left", "yanchor":"top"},
            "margin": dict(l=80, r=80, t=10, b=20),
        },
    )
    fig_hydro.update_layout(xaxis_title="Porosity []", yaxis_title="Predicted permeability [grain size input unit^2]")
    fig_hydro.update_yaxes(type="log")
    
    #plot
    data = [go.Scatter(**plot) for plot in graph_data]
    data += [go.Scatter(x=xdata, y=ydata, mode='markers', marker={"color":"blue"}, name="Data")]
    fig = go.Figure(
        data=data,
        layout = {
            "legend": {"xanchor":"right", "yanchor":"top"},
            "margin": dict(l=80, r=80, t=10, b=20),
        },
    )
    fig.update_xaxes(type="log")
    fig.update_layout(xaxis_title=header[0], yaxis_title=header[1])
    
    return fig, children_fitted_params, children_geo_params, fig_hydro


@app.callback(
    Output("download-results", "data"),
    Input("download-results-btn", "n_clicks"),
    State("dropdown-out-format-selector", "value"),
    State("fitted-result-data", "children"),
    prevent_initial_call=True,
)
def download_results(n_clicks, fmt, fit_res_table):
    #reparse res_table
    #this is awkward because the table object deeply nested...
    the_table = fit_res_table['props']['children'][0]
    column_header = [
        x['props']['children'] for x in the_table['props']['children'][0]['props']['children']['props']['children']
    ]
    #note: parameter are writted sequentially in the table, so if we pass value to the function in the order we are fine!
    fit_params = []
    for nested_dict in the_table['props']['children'][1]['props']['children']:
        params = []
        for i in range(1,len(nested_dict['props']['children'])):
            params.append(float(nested_dict['props']['children'][i]['props']['children']))
        fit_params.append(params)
    fit_params = [[np.log10(x[i]) for x in fit_params] for i in range(len(fit_params[0]))]
    #plot function
    x_th = np.logspace(fit_params[0][4],6, 1000)
    y_th = [fitter.Fredlund_PSD(x_th, *params) for params in fit_params]
    df = pd.DataFrame(
        {'Grain size': x_th, **{column_header[i+1]:y_th[i] for i in range(len(fit_params))}}
    )
    #send!
    output = io.BytesIO()
    if fmt == "CSV (.csv)":
        df.to_csv(output, index=False)
    elif fmt == "Excel (.xlsx)":
        df.to_excel(output, index=False)
    data = output.getvalue()
    return dcc.send_bytes(data, "calibrated_PSD."+fmt.split('.')[1].split(')')[0])

if __name__ == "__main__":
    app.run_server(debug=True)
