import csv

file = open("./summary_run3_fams.txt", 'r')

chart = []
gene_to_col = {}

for line, info in enumerate(file):
    temp_info = info.strip() # remove the "\n" at the end of the string

    # determine what the line of chart contains
    #   1. species
    #       - a line has a species if there's a colon AND it does not start with "<"
    #   2. node
    #       - ignore for now
    #       - a node is on a line with a colon and the line starts with "<"
    #   3. gene families
    #       - if there is a comma in the line

    quit = False

    temp_row = []
    if line is 1:
        # set up the columns (gene families)

        temp_row.append("species")
        all_genes = temp_info.split(',') # give me everything that is separated by commas

        for g in all_genes:
            temp_row.append(g)
            gene_to_col[g] = len(temp_row) - 1 # save the column the gene is stored in for later

        chart.append(temp_row)
    elif line > 1:
        # process species and nodes

        # is it a species?
        if (":" in info) and (info[0] is not "<"):
            temp_row = [0 for i in range(len(chart[0]))]   # pre-allocate the next row so that it can be filled in by counts later

            species_info = temp_info.split(':')

            species_name = species_info[0]
            species_genes = species_info[1]             # right now, species_genes is a string
            species_genes = species_genes.split(',')    # now, species genes is a list of strings
            species_genes[0] = species_genes[0][1:]     # remove the tab character, which is the first one in the first character in the first string in species_genes

            temp_row[0] = species_name[:species_name.find('<')]

            # assign the genes to their respective columns
            for g in species_genes:
                # get the family and count
                separator = g.find('[')
                family = g[:separator]
                count = g[separator+1:-2]   # we don't want the * or ]

                # figure out where the count should go in the next row based on the family
                temp_row[gene_to_col[family]] = count

            # add the next row of species and genes to the chart
            chart.append(temp_row)

# make sure the chart is the right size
print(len(chart), len(chart[0]))


# export the chart to csv
file = open('gene_family_counts.csv', 'w+', newline='')
with file:
    write = csv.writer(file)
    write.writerows(chart)
