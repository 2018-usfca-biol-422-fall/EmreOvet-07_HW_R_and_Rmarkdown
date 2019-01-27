Analysis of BLAST Results
================
Emre Ovet
October 12, 2018

Introduction
============

The human skin contains a large amount of bacteria on its surface due to its constant contact with air and different objects in daily life. Those bacterial communities are able to survive on skin surface for long periods of time because most of them are highly resistant to environmental stress factors like temperature and UV, thus, skin surface is a favorable habitat for them.

According to the recent studies, it is found out that the diversity of bacterial communities on palm surface is much higher than what is previously assumed and that after washing hands, the surviving bacteria are able to reproduce rapidly and recover their previous population. It is also found out that the genera of bacteria living on the palm surface largely varies between individuals.

In the light of these information, we hypothesize that the specificity of the bacterial community genera in male and female palm surfaces can be used for forensic investigation in terms of gender identification by comparing the bacteria found on different objects with those found on the palm surfaces of individuals. In order to test our hypothesis, we made several plots using the study of Fierer et al. in which bacterial communities from keyboards and computer mouses were extracted along with those from the palm surfaces of individuals using those equipment, then the genera of the bacteria from keyboards and computer mouses were compared with those from skin surfaces using pyrosequencing. The first figure will demonstrate the most abundant taxa in male vs female hands to see the level of similarities and differences between bacterial community genera in male and female palm surfaces, the second and third figures will demonstrate the of the abundance of the most populous bacterial taxa by gender in sebum and dust environments to see if those two bacterial taxa present equally in those surfaces, the fourth and fifth figures will demonstrate the host subject identities of the individuals with the most populous bacterial taxa by gender to see if the most populous bacterial taxa are present equally among the individuals with same gender.

Methods
=======

Sample origin and sequencing
----------------------------

For the computer mouse analysis, the palm surfaces of 5 men and 4 women from University of Colorado campus were swabbed along with the computer mouses of each individual's personal computer which were last touched by their owners 12 hours before the swabbing. All swabs were stored in -80°C until the DNA extraction. Then, the database of bacterial communities from both the right and left hand palm surfaces of 270 other individuals were compiled, sampled and analyzed.

For the pyrosequencing, DNA was extracted from samples by using MO BIO Powersoil DNA isolation kit. Amplicons were cleaned using the UltraClean-htp 96-well PCR Clean-up kit and the DNA concentration inside the amplicon were measured using the Quant-iT PicoGreen dsDNA reagent and kit.

Computational
-------------

On a tmux session, QC reports were created using the data from NCBI using FastQC program. Then, we used the Cyberduck program to transfer the files from tmux to our computer. After then, we trimed the sequences by looking at their quality scores, that is, any sequence below 150 base pairs were filtered. Furthermore, we Blasted the sequences to obtain the top match of each sequence.Finally, we used dplyr functions to visualize our data on Rstudio by making plots.

Results
=======

We found out that there are only 5 bacterial taxa that are widely common among males and just 1 in females which is not expected. *Solemya pervernicosa gill symbiont* is the most populous bacterial taxa among the 5 most abundant in males and *Bartonella washoensis* is the most populous among the females. Both of these bacterial taxa are abundant on both sebum and dust surfaces ,with *Solemya pervernicosa gill symbiont* being more abundant in dust and *Bartonella washoensis* in sebum. Moreover, we found out that 1 female individual out of 4 possessed all the *Bartonella washoensis* on her which was not expected. We obtained a higher amount of male individuals possessing *Solemya pervernicosa gill symbiont* with 3 males out of 5. Overall, we can say that our hypothesis is true.

``` r
# Be sure to install these packages before running this script
# They can be installed either with the install.packages() function
# or with the 'Packages' pane in RStudio

# load packages
library("dplyr")
library("tidyr")
library("knitr")
library("ggplot2")
```

``` r
# Output format from BLAST is as detailed on:
# https://www.ncbi.nlm.nih.gov/books/NBK279675/
# In this case, we used: '10 sscinames std'
# 10 means csv format
# sscinames means unique Subject Scientific Name(s), separated by a ';'
# std means the standard set of result columns, which are:
# 'qseqid sseqid pident length mismatch
# gapopen qstart qend sstart send evalue bitscore',


# this function takes as input a quoted path to a BLAST result file
# and produces as output a dataframe with proper column headers
# and the 'qseqid' column split into sample and seq number
read_blast_output <- function(filename) {
  data_in <- read.csv(filename,
                      header = FALSE, # files don't have column names in them
                      col.names = c("sscinames", # unique Subject Sci Name(s)
                                    "qseqid",    # Query Seq-id
                                    "sseqid",    # Subject Seq-id
                                    "pident",    # Percntge of identical matches
                                    "length",    # Alignment length
                                    "mismatch",  # Number of mismatches
                                    "gapopen",   # Number of gap openings
                                    "qstart",    # Start of alignment in query
                                    "qend",      # End of alignment in query
                                    "sstart",    # Start of alignment in subj
                                    "send",      # End of alignment in subject
                                    "evalue",    # Expect value
                                    "bitscore"))  # Bit score

  # Next we want to split the query sequence ID into
  # Sample and Number components so we can group by sample
  # They originally look like "ERR1942280.1"
  # and we want to split that into two columns: "ERR1942280" and "1"
  # we can use the separate() function from the tidyr library to do this
  # Note that we have to double escape the period for this to work
  # the syntax is
  # separate(column_to_separate,
  # c("New_column_name_1", "New_column_name_2"),
  # "seperator")
  data_in <- data_in %>%
    separate(qseqid, c("sample_name", "sample_number"), "\\.")
}
```

``` r
# this makes a vector of all the BLAST output file names, including
# the name(s) of the directories they are in
files_to_read_in <- list.files(path = "output/blast",
                               full.names = TRUE)

# We need to create an empty matrix with the right number of columns
# so that we can rbind() each dataset on to it
joined_blast_data <- matrix(nrow = 0,
                            ncol = 14)

# now we loop over each of the files in the list and append them
# to the bottom of the 'joined_blast_data' object
# we do this with the rbind() function and the function we
# made earlier to read in the files, read_blast_output()
for (filename in files_to_read_in) {
  joined_blast_data <- rbind(joined_blast_data,
                             read_blast_output(filename))
}
```

``` r
# Next we want to read in the metadata file so we can add that in too
# This is not a csv file, so we have to use a slightly different syntax
# here the `sep = "\t"` tells the function that the data are tab-delimited
# and the `stringsAsFactors = FALSE` tells it not to assume that things are
# categorical variables
metadata_in <- read.table(paste0("data/metadata/",
                                 "fierer_forensic_hand_mouse_SraRunTable.txt"),
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)

# Finally we use the left_join() function from dplyr to merge or 'join' the
# combined data and metadata into one big table, so it's easier to work with
# in R the `by = c("Run_s" = "sample_name")` syntax tells R which columns
# to match up when joining the datasets together
joined_blast_data_metadata <- metadata_in %>%
  left_join(joined_blast_data,
            by = c("Run_s" = "sample_name"))
```

``` r
joined_blast_data_metadata %>%
  group_by(sscinames, sex_s) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(sex_s != "Not applicable") %>%
  filter(n > 100) %>%
  ggplot(aes(sscinames, y = n, fill = sex_s)) +
  geom_col() +
  theme(axis.text = element_text(angle = 90,
                                 hjust = 1)) +
  ggtitle("Number of most abundant taxa in male vs female hands")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github/Different%20bacterial%20taxa%20on%20male%20vs%20female%20hands-1.png)

**Figure 1** In this figure, we can see that the number of most abundant taxa in male and female hands which are higher than 100. There are only 7 bacterial taxa out of hundreds that are widely common, and only 2 of them are found in female palm surfaces; Acidovorax sp. and Bartonella washoensis. Overall, Solemya pervemicosa gill symbiont is the most abundant bacteria in males and Bartonella washoensis is the most abundant bacteria in females.

``` r
joined_blast_data_metadata %>%
  group_by(sscinames, env_material_s) %>%
  tally() %>%
  filter(sscinames %in% c("Bartonella washoensis")) %>%
  ggplot(aes(sscinames, y = n, fill = env_material_s)) +
  geom_col(position = position_dodge()) +
    theme(axis.text = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundance of Bartonella washoensis in male vs female hands")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github/Abundance%20of%20Bartonella%20washoensis%20in%20dust%20vs%20sebum-1.png)

**figure 2** In this figure, we can see that the abundance of *Bartonella washoensis* is between 750 and 900 in dust and between 500 and 750 in sebum.

``` r
joined_blast_data_metadata %>%
  group_by(sscinames, env_material_s) %>%
  tally() %>%
  filter(sscinames %in% c("Solemya pervernicosa gill symbiont")) %>%
  ggplot(aes(sscinames, y = n, fill = env_material_s)) +
  geom_col(position = position_dodge()) +
    theme(axis.text = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundance of Solemya pervernicosa gill symbiont in dust vs sebum")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github/Abundance%20of%20Solemya%20pervernicosa%20gill%20symbiont%20in%20dust%20vs%20sebum-1.png)

**figure 3** In this figure, we can see that the abundance of *Solemya pervernicosa gill symbiont* is between 750 and 1000 in dust and between 1500 and 1750 in sebum.

``` r
joined_blast_data_metadata %>%
  group_by(sscinames, host_subject_id_s) %>%
  tally() %>%
  filter(sscinames %in% c("Solemya pervernicosa gill symbiont")) %>%
  ggplot(aes(sscinames, y = n, fill = host_subject_id_s)) +
  geom_col(position = position_dodge()) +
    theme(axis.text = element_text(angle = 90, hjust = 1)) +
  ggtitle("Host subject identities of the
          individuals with Solemya pervernicosa gill symbiont")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github/Host%20subject%20identities%20of%20the%20individuals%20with%20Solemya%20pervernicosa-1.png)

**figure 4** In this figure, we can see that there are many individuals that have *Solemya pervernicosa gill symbiont* on their palm surfaces. The indivual with host identity M2 seems to be a female since the abundance of the *Solemya pervernicosa gill symbiont* on their palms are close to zero. It is not certain that the rest of the individuals are males since we know that there are individuals with unidentified genders that carry the *Solemya pervernicosa gill symbiont*. In order to figure out how many of the males possess the *Solemya pervernicosa gill symbiont*, we have to compare this graph with another graph showing host subject identities of the individuals with *Bartonella washoensis*.

``` r
joined_blast_data_metadata %>%
  group_by(sscinames, host_subject_id_s) %>%
  tally() %>%
  filter(sscinames %in% c("Bartonella washoensis")) %>%
  ggplot(aes(sscinames, y = n, fill = host_subject_id_s)) +
  geom_col(position = position_dodge()) +
    theme(axis.text = element_text(angle = 90, hjust = 1)) +
  ggtitle("Host subject identities of the
          individuals with Bartonella washoensis")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github/Host%20subject%20identities%20of%20the%20individuals%20with%20Bartonella%20washoensis-1.png)

**figure 5**

In this figure, we can see that the individuals with host identities F6, F7 and M8 are the ones with unidentified genders since they both show up in both graphs with high abundances. We can also identify the individuals F2, F8 and M7 as males since their abundances of *Bartonella washoensis* are close to zero.

``` r
levels(factor(joined_blast_data_metadata$anonymized_name_s))
```

    ##  [1] "F2Mose396" "F2Plmr396" "F5Mose396" "F5Plmr396" "F6Mose396"
    ##  [6] "F6Plmr396" "F7Mose396" "F7Plmr396" "F8Mose396" "F8Plmr396"
    ## [11] "M1Mose396" "M1Plmr396" "M2Mose396" "M2Plmr396" "M7Mose396"
    ## [16] "M7Plmr396" "M8Mose396" "M8Plmr396" "M9Mose396" "M9Plmr396"

``` r
joined_blast_data_metadata %>%
  group_by(anonymized_name_s) %>%
  summarize(mean_pident = mean(pident),
            sd_pident = sd(pident)) %>%
  ggplot(aes(x = anonymized_name_s,
             y = mean_pident)) +
  geom_col(fill = "dodger blue") +
  geom_errorbar(aes(ymax = mean_pident + sd_pident,
                    ymin = mean_pident - sd_pident),
                width = 0.3) +
theme(axis.text.x = element_text(angle = 90,
                                    hjust = 1))
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github/Mean%20percent%20identity%20of%20query%20sequences-1.png)

``` r
# Here we're using the dplyr piping syntax to select a subset of rows matching a
# criteria we specify (using the filter) function, and then pull out a column
# from the data to make a histogram.
joined_blast_data_metadata %>%
  filter(env_material_s == "dust") %>%
  filter(grepl("F2", host_subject_id_s)) %>%
  ggplot(aes(x = pident)) +
    geom_histogram() +
    ggtitle("Percent Identity") +
    xlab("Percent")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Analysis_of_BLAST_Results_files/figure-markdown_github/histograms-1.png)

``` r
# Finally, we'd like to be able to make a summary table of the counts of
# sequences for each subject for both sample types. To do that we can use the
# table() function. We add the kable() function as well (from the tidyr package)
# in order to format the table nicely when the document is knitted
kable(table(joined_blast_data_metadata$host_subject_id_s,
            joined_blast_data_metadata$sample_type_s))
```

|     |  computer mouse|  right palm|
|-----|---------------:|-----------:|
| F2  |             396|         410|
| F5  |             365|         777|
| F6  |             662|         422|
| F7  |             655|         546|
| F8  |             878|         374|
| M1  |             456|         878|
| M2  |             670|         775|
| M7  |             970|         689|
| M8  |             717|         280|
| M9  |             571|         968|

Discussion
==========

We found out that Solemya pervemicosa gill symbiont is the most abundant bacteria in males and Bartonella washoensis is the most abundant bacteria in females. We're wondering why this is the case but it is hard to find an answer since the information about those bacteria are very limited for now. As far as we know, Solemya pervemicosa gill symbionts are chemosynthetic bacteria, that is, they convert CH4 or CO2 into organic compounds and they provide their host with necessary nutrients. We can argue that these symbionts were somehow transferred from Solemya pervemicosa clams to humans and they had a very favorable habitat in palm surface, specifically male palm surfaces since none of them were found in female palm surfaces. One possibility is that a male individual made contact with the gills of a Solemya pervemicosa and then he made its way to the men's restroom in university campus and from here, it has spread to the hands of the male individuals from siphon and faucet surfaces. Further experimentation and research is needed to answer how did they transfer to humans and why they are only seen in males in a very high amount.

When it comes to Bartonella washoensis, we have to urge upon a detail, that is, Bartonella washoensis is pathogenic for humans as it causes meningitis. Somehow, this pathogenic bacteria didn't make its way into the body and cause disease because we know that all the female individuals participated in this experiment are healthy. One possibility is that they are all immunized by meningococcal vaccine. Another interesting detail in our results is that this bacteria is only found in female palm surfaces. We know that this bacteria is transmitted to humans by ground squirrels, so we can argue that the female individual with host identity M2 made contact with a ground squirrel, leading the displacement of some of that bacteria from squirrel's fur to her palm surface and the bacteria were able to survive and reproduce. Once again, further experimentation and research is needed to answer why those female individuals didn't have meningitis and why they are only seen in females in a very high concentration.

If we can find out the answers of those questions and detect other bacteria that are gender-specific, it'll be useful in forensic investigations because we can find out if the suspect is a male or female by analyzing the bacteria genera and determining the gender-specific ones.
