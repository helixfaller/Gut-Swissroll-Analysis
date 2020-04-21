#!/usr/bin/env python
# coding: utf-8

# In[7]:

def box_plot(columns, genotype, data, test="Mann-Whitney",
             hue="Sacrifice Age (month)",hue_split=0, loc="inside",
            genotype_column="Genotype", palette="Set2_r"):
    
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt
    import seaborn as sns
    from statannot import add_stat_annotation
    from itertools import combinations
    plt.rcParams["figure.figsize"] = 12.0, 9.0
    plt.rcParams["xtick.labelsize"] = 15
    plt.rcParams["ytick.labelsize"] = 15
    
    ## get boxpairs depending on hue and hue_split
    liste_box_pairs = []
    if hue == None:
        for i in range(len(genotype)):
            variable = (genotype[i])
            liste_box_pairs.append(variable)
    elif hue != None:
        splitting = data[hue].tolist()
        splitting_variables = []
        if hue_split != None:
            for variable in splitting:
                if variable not in splitting_variables:
                    splitting_variables.append(variable)
                else:
                    pass
            for x in range(len(genotype)):
                variable = (genotype[x], splitting_variables[hue_split])
                liste_box_pairs.append(variable) 
        elif hue_split == None:
            for variable in splitting:
                if variable not in splitting_variables:
                    splitting_variables.append(variable)
                else:
                    pass
            for i in range(len(splitting_variables)):
                for x in range(len(genotype)):
                    variables = (genotype[x], splitting_variables[i])
                    liste_box_pairs.append(variables)
        else:
            pass
    else:
        pass
    box_pairs = list(combinations(liste_box_pairs, 2))

    ax = sns.boxplot(data=data, x=genotype_column, y=columns, hue=hue,
                         linewidth=4, palette=palette, fliersize=6)
    add_stat_annotation(ax, data=data,x=genotype_column, y=columns,
                            hue=hue, box_pairs=box_pairs, text_format="star",
                            test=test, loc=loc, verbose=5)
    ax.set(xlabel=None)
    ax.set_ylabel(columns, fontsize=15)
    # if max(self.data[self.columns].tolist()) > 500:
            #ax.set_yscale('log')
    ax.legend(title=hue, title_fontsize="x-large", fontsize="x-large",
              loc='center right', bbox_to_anchor=(1.3, 0.5), ncol=1)
    return plt.show()


def weight_curve(data, hue="Genotype", palette="Set2_r"):
    
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt
    import seaborn as sns
    from statannot import add_stat_annotation
    from itertools import combinations
    
    plt.rcParams["figure.figsize"]=7.0, 9.0
    
    ## get weight columns (columns where first four letters are "Week"
    weight_columns = []
    for columns in data.columns.tolist():
        if columns[:4] == "Week":
            weight_columns.append(columns)
            
    ## get for each mouse, timpoint and weight of this thimepoint
    mice = []
    times = []
    time_weights = []
    hues = []
    for mouse in data.index.tolist():
        for a, time in enumerate(weight_columns):
            time_weight = data[weight_columns].loc[mouse][a]
            mice.append(mouse)
            times.append(time)
            time_weights.append(time_weight)
            hue_ = data[hue][mouse]
            hues.append(hue_)
    
    ## new dataframe
    data_df2 = {"Mouse": [], "Time": [], "Time_weight": []}
    df2 = pd.DataFrame(data_df2)
    
    df2["Mouse"] = mice
    df2["Time"] = times
    df2["Time_weight"] = time_weights
    df2[hue] = hues
    
    ## replace weeks with number i
    for i, column in enumerate(weight_columns):
        df2 = df2.replace(column, i)
    
    ## get nice graph
    ## get nice graph
    ax = sns.pointplot(x="Time", y="Time_weight", data=df2, capsize=.2,
                      hue=hue,eight=6, aspect=.75,
                      palette=palette, ci=None)
    ax.set_xlabel('time [weeks]', fontsize = 15)
    ax.set_ylabel('weight [g]', fontsize = 15)
    return plt.figure()

def weightloss(data, time):
    
    import pandas as pd
    import seaborn as sns
    from matplotlib import pyplot as plt
    plt.rcParams["figure.figsize"] = 12.0, 9.0
    
    ## get weight columns, that are selected by the string time
    for column in data.columns.tolist():
        if column[:4] != time:
            data = data.drop([column], axis="columns")
        else:
            pass
    
    ## fill all NaN values with 0
    data = data.fillna(0)
    
    ## get max weight and last weight of each mouse
    max_weight = []
    last_weight = []
    mice = []
    for i, mouse in enumerate (data.index.tolist()):
        max_weight.append(max(data.loc[mouse]))
        last_weight.append(data.loc[mouse][len(data.columns.tolist())-1])
        mice.append(mouse)

    ## calculate the percentage weight loss
    percentage_weight_loss = [(max_weight - last_weight)/(max_weight) for max_weight,
                              last_weight in zip(max_weight, last_weight)]
    
    ## new dataframe with max weight and last weight; index is mouse ID
    data_df2 = {"max weight": [], "last weight": [], "percentage weight loss": []}
    df2 = pd.DataFrame(data_df2)
    df2["max weight"] = max_weight
    df2["last weight"] = last_weight
    df2["percentage weight loss"] = percentage_weight_loss
    df2.index = data.index.tolist()
    
    ## if last weight == 0; drop row, because mouse was already euthanized
    for mouse in df2.index.tolist():
        if df2.loc[mouse]["last weight"] == 0:
            df2 = df2.drop(index=mouse)
        else:
            pass

    ## drop all mice, where percentage weight loss is below 10%, since mice become euthanized
    ## after a weight loss of 20%
    for mouse in df2.index.tolist():
        if df2.loc[mouse]["percentage weight loss"] < 0.1:
            df2 = df2.drop(index=mouse)
        else:
            pass
            

    ## get nice scatterplot
    ## hue is a grouping variable that will produce points with different colors
    ## s represents the radius of the point
    ax = sns.scatterplot(x="max weight", y="percentage weight loss", hue=df2.index, data=df2, s=100)
    ## set y min to 0.10 == 10%; data <10% is irrelevant
    ax.set_ylim(ymin=0.10)
    ## set y label
    ax.set_ylabel("weight loss [%]", fontsize=15)
    ## set x label
    ax.set_xlabel("max weight [g]", fontsize=15)
    ## add legend
    ax.legend(title_fontsize="x-large", fontsize="x-large",
              loc='center right', bbox_to_anchor=(1.3, 0.5), ncol=1)
    return plt.figure()