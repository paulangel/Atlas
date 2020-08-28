import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
print(dcc.__version__)
import socket
socket.gethostbyname('localhost')

import numpy, scipy, sys, os, pandas, sklearn, scanpy, sklearn.decomposition, sklearn.cluster, gc
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *
import functions

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

output_loc  = '/Users/pwangel/Downloads/'
data        = pandas.read_csv('/Users/pwangel/Atlas/data/blood_atlas_expression.tsv', sep='\t', index_col=0)
annotations = pandas.read_csv('/Users/pwangel/Atlas/data/blood_atlas_annotations.tsv', sep='\t', index_col=0) 

data = functions.transform_to_percentile(data)
genes = functions.calculate_platform_dependence(data, annotations)

blood_atlas_colours = pandas.read_csv('/Users/pwangel/Atlas/data/blood_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}

available_indicators = ['Platform_Category', 'Dataset', 'tier1', 'celltype']

app.layout = html.Div([
    html.Div([

        html.Br(),
        html.Label('Atlas Data'),
        html.Div([
            dcc.Dropdown(
                id='Population-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='tier1'
            ),
        ],
        style={'width': '100%'}), html.Br(),
    ], className='three columns'),

    html.Div([
    dcc.Graph(id='atlas-graphic')
    ], 
    className='nine columns'),

], className='row')

#app.css.append_css({
#    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
#})

@app.callback(
    Output('atlas-graphic', 'figure'),
    [Input('Population-column', 'value')
     ], 
     state=[dash.dependencies.State('atlas-graphic', 'figure')])
def update_pca(population_column, previous_figure):

    df_anno = annotations
    df_data = data[df_anno.index] 

    pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
    pca_coords = pca.fit_transform(df_data.loc[genes.Platform_VarFraction.values<=0.2].transpose())

    fig = Figure()
   
    if df_anno[population_column].dtype==numpy.float64():
        fig.add_trace(Scatter3d(x=pca_coords[:,0], y=pca_coords[:,1], z=pca_coords[:,2], 
                    mode='markers', text=df_anno.display_metadata.values, 
                    opacity=0.9, name='Population',
                    marker=dict(size=5, color=df_anno[population_column].values,
                    colorscale = 'RdBu', colorbar = dict(x=0.85))))
    else:
        for i_type in df_anno[population_column].unique():
            if i_type in blood_atlas_colours.keys():
                i_colour = blood_atlas_colours[i_type]
            else:
                i_colour = numpy.random.choice(['red', 'blue', 'green', 'black'])
                print(i_type)
            sel = df_anno[population_column].values == i_type 
            fig.add_trace(Scatter3d(x=pca_coords[sel,0], y=pca_coords[sel,1], z=pca_coords[sel,2], 
                    mode='markers', text=df_anno.display_metadata.values[sel], 
                    opacity=0.9, name=i_type, 
                    marker=dict(size=5, color=i_colour)))

    fig.update_layout(uirevision=True,width=1000,height=630,
        scene = dict(xaxis_title='PC1 %%%.2f' %(pca.explained_variance_ratio_[0]*100),
        	       yaxis_title='PC2 %%%.2f' %(pca.explained_variance_ratio_[1]*100),
        	       zaxis_title='PC3 %%%.2f' %(pca.explained_variance_ratio_[2]*100))
        )

    plot(fig, auto_open=False, filename=output_loc+'/atlas.html')

    return fig

if __name__ == '__main__':
    #app.run_server(mode='inline', host = '127.0.0.1')
    app.run_server(host = '127.0.0.1')




