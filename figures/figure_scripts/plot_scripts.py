####
import altair as alt
import pandas as pd

###POWER CURVE - VARIANT LEVEL ###
def create_df(highest_rank, compare, successTable):
    x = range(1,highest_rank)
    denom = len(successTable)
    data=[]
    for runtype in compare:
        success = successTable.loc[(successTable['Variant_Level_noMOI_'+str(runtype)]=='Variant_Present_noMOI') | (successTable['Variant_Level_noMOI_'+str(runtype)]=='Variant_Present_noMOI')]
        for i in x:
            num = len(success.loc[success['Variant_Level_noMOI_rank_'+str(runtype)] <= i])
            data.append([i,num, (num/denom)*100, runtype])
        try:
            print(runtype,'highest rank:', max(success['Variant_Level_noMOI_rank_'+str(runtype)]), 'max%:', (num/denom)*100, num)
        except:
            print(runtype, '0 percent success')

    df = pd.DataFrame(data, columns=['Rank', 'NumPatients', 'Percent_Variants', 'Run_Type'])
    return df, denom


def create_plot(denom, source, compare, color_scheme, big_plot,domain, miser='Exomiser'):
    title=str(denom) +'Diagnostic Variants (Variant Level - no MOI requirement)'
    if big_plot:
        bigChart = alt.Chart(source, title=title).mark_line().encode(
            x=alt.X('Rank', title=str(miser)+' Rank of Causal Variant'),
            y=alt.Y('Percent_Variants', title='Percent of Causal Variants within ' +str(miser)+' Rank', scale=alt.Scale(domain=[0,100])),
            color=alt.Color('Run_Type:N',sort=compare,scale=alt.Scale(scheme=color_scheme, domain=domain), legend=alt.Legend(orient='top')),
            tooltip=['Rank', 'Percent_Variants', 'Run_Type']
        )
    
    
    zoom_source = source.loc[source['Rank'] <=30]

    zoomChart = alt.Chart(zoom_source, title=title).mark_line(point=alt.OverlayMarkDef(size=50)).encode(
        x=alt.X('Rank', title=str(miser)+' Rank of Causal Variant'),
        y=alt.Y('Percent_Variants', title='Percent of Causal Variants within '+str(miser)+' Rank', scale=alt.Scale(domain=[0,100])),
        color=alt.Color('Run_Type:N', sort=compare, scale=alt.Scale(scheme=color_scheme, domain=domain)),
        tooltip=['Rank', 'Percent_Variants', 'Run_Type']
    )
    if not big_plot:
        plot = zoomChart
    else:
        print('concatinating charts')
        plot=alt.hconcat(bigChart, zoomChart)#.configure_legend(labelLimit=0).configure_axis(
        # labelFontSize=15,
        # titleFontSize=15).configure_legend(labelLimit=0,labelFontSize=15, titleFontSize=15)#.configure_axis(grid=False)
    return plot


###POWER CURVE - GENE LEVEL ###
def create_df_gene(highest_rank, compare, successTable):
    x = range(1,highest_rank)
    denom = len(successTable)
    data=[]
    for runtype in compare:
        success = successTable.loc[(successTable['Gene_Level_'+str(runtype)]=='Exomiser_Success')]
        for i in x:
            num = len(success.loc[success['Gene_Level_rank_'+str(runtype)] <= i])
            data.append([i,num, (num/denom)*100, runtype])
        print(runtype, 'highest rank:', max(success['Gene_Level_rank_'+str(runtype)]), 'max%:', (num/denom)*100, num)

    df = pd.DataFrame(data, columns=['Rank', 'NumPatients', 'Percent_Variants', 'Run_Type'])
    return df, denom


def create_plot_gene(denom, source, compare, color_scheme, big_plot, domain,miser='Exomiser'):

    title=str(denom) +'Diagnostic Genes (Gene Level)'
    if big_plot:
        bigChart = alt.Chart(source, title=title).mark_line().encode(
            x=alt.X('Rank', title=str(miser)+' Rank of Causal Gene'),
            y=alt.Y('Percent_Variants', title='Percent of Causal Genes within ' +str(miser)+ ' Rank', scale=alt.Scale(domain=[0,100])),
            color=alt.Color('Run_Type:N',sort=compare,scale=alt.Scale(scheme=color_scheme)),
            tooltip=['Rank', 'Percent_Variants', 'Run_Type']
        ).properties(
            width=600,
            height=500)
    
    
    zoom_source = source.loc[source['Rank'] <=30]

    zoomChart = alt.Chart(zoom_source, title=title).mark_line(point=alt.OverlayMarkDef(size=50)).encode(
        x=alt.X('Rank', title=str(miser)+' Rank of Causal Gene'),
        y=alt.Y('Percent_Variants', title='Percent of Causal Genes within ' +str(miser)+ ' Rank', scale=alt.Scale(domain=[0,100])),
        color=alt.Color('Run_Type:N', sort=compare, scale=alt.Scale(scheme=color_scheme, domain=domain)),
        tooltip=['Rank', 'Percent_Variants', 'Run_Type']
    ).properties(
        width=600,
        height=500)
    
    

    if not big_plot:
        plot = zoomChart
    else:
        plot=alt.hconcat(bigChart, zoomChart)#.configure_legend(labelLimit=0).configure_axis(
        # labelFontSize=15,
        # titleFontSize=15).configure_legend(labelLimit=0,labelFontSize=15, titleFontSize=15)#.configure_axis(grid=False)
    return plot

order=['NONE', 'ALPHA_MISSENSE', 'MVP', 'REVEL', 'SPLICE_AI', 'CADD','MUTATION_TASTER', 'POLYPHEN', 'SIFT','REMM', 'Not_Prioritized']
def create_score_plot(run_type, table, order):
    base = alt.Chart(table, title=str(len(table))+' Variants from UDN Individuals(' + str(run_type)+')')
    bars=base.mark_bar().encode(
        y=alt.Y('variant_type:N', sort='-x', title=None).axis(offset=5, domainOpacity=0),
        x=alt.X('count()', title='Number of Diagnostic Variants'), 
        color=alt.Color('Variant_Level_noMOI_maxPathSource_' + str(run_type), sort=order,scale=alt.Scale(domain=order,scheme='tableau20')),
        tooltip = ['count()','variant_type:N','Variant_Level_noMOI_maxPathSource_' + str(run_type),],
    )

    text = base.mark_text(
        align='left',
        baseline='middle',
        dx=3,
        color='black',
        size=12
    ).encode(
        y=alt.Y('variant_type:N', sort='-x', title=None).axis(offset=5, domainOpacity=0),
        x=alt.X('count()'), 
        text='count():Q'
    )
    plot = alt.layer(bars, text).resolve_scale(color='independent').properties(height=150, width=300)
    return plot

##DOT PLOT##
def dot_df_pvals(source, compare):
    data = []
    p_values = [1, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01]
    for runType in compare:
        print(runType)
        table = source[['Variant_Level_noMOI_'+runType,'Variant_Level_noMOI_p_'+runType]]

        successes = table.loc[table['Variant_Level_noMOI_'+runType]=='Variant_Present_noMOI']
        all_fails = len(table.loc[table['Variant_Level_noMOI_'+runType].isin(['Variant_Not_Present_noMOI', 'Gene_Not_in_Output'])])
        all_successes = len(successes)
        denom = len(table)
        print(all_successes + all_fails, 'total genes')
        for p in p_values:
            success = len(successes.loc[successes['Variant_Level_noMOI_p_'+runType]<=p])
            fail = len(successes.loc[(~successes['Variant_Level_noMOI_p_'+runType].notnull()) | (successes['Variant_Level_noMOI_p_'+runType] > p)])
            data.append([runType, (success/denom)*100, p])
        data.append([runType, (all_fails/denom)*100, 'failed'])
    df = pd.DataFrame(data, columns=['RunType', 'SuccessRate', 'SuccessType'])
    plot = alt.Chart(df).mark_circle(size=40, filled=False, opacity=1).encode(
        x=alt.X('SuccessRate', sort=compare,title='Percent of True Causal Variants Within Exomiser P' ,scale=alt.Scale(domain=[0,100])),
        y= alt.Y('RunType', ),
        color=alt.Color('SuccessType:N', sort='-x', scale=alt.Scale(scheme='category10')),
        tooltip = ['RunType', 'SuccessRate', 'SuccessType']
    ).properties(
        height=200,
        width=600
    ).configure_legend(labelLimit=0).configure_axis(
        labelFontSize=15,
        titleFontSize=15).configure_scale(
        bandPaddingInner=0.9

    )
    return plot

def dot_df_ranks(source, compare, color_scheme='category10'):
    data = []
    ranks = ['All', 30, 20, 10, 5, 1]
    for runType in compare:
        print(runType)
        table = source[['Variant_Level_noMOI_'+runType,'Variant_Level_noMOI_rank_'+runType]]

        successes = table.loc[table['Variant_Level_noMOI_'+runType]=='Variant_Present_noMOI']
        all_fails = len(table.loc[table['Variant_Level_noMOI_'+runType].isin(['Variant_Not_Present_noMOI', 'Gene_Not_in_Output'])])
        all_successes = len(successes)
        denom = len(table)
        print(all_successes + all_fails, 'total variants')
        for rank in ranks:
            if rank=='All':
                data.append([runType, (all_successes/denom)*100, 'Success'])
            else:
                success = len(successes.loc[successes['Variant_Level_noMOI_rank_'+runType]<=rank])
                data.append([runType, (success/denom)*100, rank])
        data.append([runType, (all_fails/denom)*100, 'Fail'])
    df = pd.DataFrame(data, columns=['RunType', 'SuccessRate', 'SuccessType'])
    plot = alt.Chart(df).mark_circle(size=100, filled=True, stroke='black', opacity=1).encode(
        x=alt.X('SuccessRate', sort=compare,title='Percent of Causal Variants Within Exomiser Rank' ,scale=alt.Scale(domain=[0,100])),
        y= alt.Y('RunType', sort=compare),
        color=alt.Color('SuccessType:N', sort='-y', scale=alt.Scale(scheme=color_scheme)),
        tooltip = ['RunType', 'SuccessRate', 'SuccessType']
    ).properties(
        height=500,
        width=100)
    # ).configure_legend(labelLimit=0).configure_axis(
    #     labelFontSize=15,
    #     titleFontSize=15).configure_scale(
    #     bandPaddingInner=0.9

    # )
    return plot


