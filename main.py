import base64
import datetime
import io

from dash import Dash, html, dcc, callback, Output, Input, State
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import flask

import fitter


instructions = [
    "1. Prepare a CSV file (comma separated) the first column being the sieve (in millimeter) and the second the percentage of passing particle (between 0 and 100)",
    "2. Drag and drop or upload the previous file. The experimental particule size distribution should be plotted on the below graph.",
    "3. Click on the button optimize to search for the fit the Fredlund model (best means with the lowest Root Mean Squared Error - RMSE).",
    "4. Fredlund model parameters, D-values and coefficient of uniformity C_U and curvature C_C are returned!",
]

layout = [
    html.H1('Particule Size Distribution Analyser', style={'textAlign':'center'}),
    html.H2('Instructions', style={'textAlign':'center'}),
] + [
    html.P(instruction, style={'textAlign':'center'}) for instruction in instructions
] + [
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select a file')
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
    html.Hr(),  # horizontal line
    html.H2('Fitted PSD curve', style={'textAlign':'center'}),
    dcc.Graph(id='graph-content'),
    html.Button('Optimize', id='optimize_button', n_clicks=0),
    html.Div(id='result-data'),
    html.Hr(),  # horizontal line
]



#TODO: download the csv and put it on /tmp
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


server = flask.Flask(__name__)
app = Dash(__name__, title="PSD-Tools", server=server)
app.layout = html.Div(layout)

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
    return fig


@app.callback(
    Output('graph-content', 'figure'),
    Output('result-data', 'children'),
    Input('optimize_button', 'n_clicks'),
    State('upload-data', 'contents'),
    State('upload-data', 'filename'),
)
def optimize(btn, contents, filename):
    if contents is None: return go.Figure(), html.Div()
    print('optimize')
    data = parse_contents(contents, filename)
    header = data.columns.tolist()
    xdata = np.array(data.iloc[:,0])
    ydata = np.array(data.iloc[:,1])
    res = fitter.fit(xdata, ydata)
    
    #plot
    x_th = np.logspace(np.log10(np.min(xdata)), np.log10(np.max(xdata)), 100)
    x_th = np.append(x_th, xdata)
    x_th.sort()
    y_th = fitter.Fredlund_PSD(x_th, *res.x)
    R2 = fitter.R2(res.fun, ydata)
    fig = go.Figure(
        data=[
            go.Scatter(x=xdata, y=ydata, mode='markers', name="Data"),
            go.Line(x=x_th, y=y_th, name="Fitted Fredlund"),
        ]
    )
    fig.update_xaxes(type="log")
    
    #compute D_10, D_60, Cu and Cc
    D_10 = fitter.Dx(10,fitter.Fredlund_PSD, res.x)
    D_30 = fitter.Dx(30,fitter.Fredlund_PSD, res.x)
    D_60 = fitter.Dx(60,fitter.Fredlund_PSD, res.x)
    Cu = D_60/D_10
    Cc = D_30**2/(D_60*D_10)
    
    #print results
    children = html.Div([
        html.P(f"RMSE (the lower the better): {np.sqrt(res.fun)}", style={'textAlign':'center'}),
        html.P(f"Determination coefficient (R2): {R2}", style={'textAlign':'center'}),
        html.P(f"a_gr: {10**res.x[0]:.3e}, n_gr: {10**res.x[1]:.3e}, m_gr: {10**res.x[2]:.3e}, dr: {10**(-res.x[3]):.6e} mm, dm: {10**(-res.x[4]):.6e} mm", style={'textAlign':'center'}),
        html.P(f"D_10: {D_10:.5f} mm, D_30: {D_30:.5f} mm, D_60: {D_60:.5f} mm", style={'textAlign':'center'}),
        html.P(f"Cu: {Cu:.1f}, Cc: {Cc:.1f}", style={'textAlign':'center'}),
    ])
    
    #predict saturated hydraulic conductivity
    
    return fig, children

if __name__ == '__main__':
    app.run_server(debug=True)
