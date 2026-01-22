#!/bin/bash

RESULTS_FILE="./gamma.out"
NUQCLQ_PATH="./edqc"
TIME_LIMIT=60

declare -A allowed_graphs
allowed_graphs=(
    ["soc-dolphins"]=1 ["rt-retweet"]=1 ["ca-netscience"]=1 ["web-polblogs"]=1
    ["rt-twitter-copen"]=1 ["email-Eu-core"]=1 ["ca-CSphd"]=1 ["web-edu"]=1
    ["ego-facebook"]=1 ["ca-GrQc"]=1 ["web-spam"]=1 ["ca-Erdos992"]=1
    ["socfb-CMU"]=1 ["ech-WHOIS"]=1 ["ca-HepPh"]=1 ["web-indochina-2004"]=1
    ["web-BerkStan"]=1 ["socfb-UCSB37"]=1 ["web-webbase-2001"]=1 ["socfb-UConn"]=1
    ["ca-AstroPh"]=1 ["socfb-UCLA"]=1 ["ca-CondMat"]=1 ["socfb-Berkeley13"]=1
    ["socfb-Wisconsin87"]=1 ["soc-epinions"]=1 ["Cit-HepTh"]=1 ["socfb-Indiana"]=1
    ["socfb-UIllinois"]=1 ["Cit-HepPh"]=1 ["socfb-UF"]=1 ["Email-Enron"]=1
    ["socfb-Penn94"]=1 ["soc-brightkite"]=1 ["socfb-OR"]=1 ["soc-slashdot"]=1
    ["wordnet-words"]=1 ["witter_combined"]=1 ["G_n_pin_pout"]=1 ["yahoo-msg"]=1
    ["web-sk-2005"]=1 ["web-uk-2005"]=1 ["rgg_n_2_17_s0"]=1 ["soc-douban"]=1
    ["wave"]=1 ["web-arabic-2005"]=1 ["Loc-Gowalla"]=1 ["soc-gowalla"]=1
    ["ca-dblp-2010"]=1 ["ca-citeseer"]=1 ["coAuthorsCiteseer"]=1 ["web-Stanford"]=1
    ["coAuthorsDBLP"]=1 ["ca-dblp-2012"]=1 ["com-dblp"]=1 ["cnr-2000"]=1
    ["web-NotreDame"]=1 ["ca-MathSciNet"]=1 ["coPapersDBLP"]=1 ["auto"]=1
    ["web-it-2004"]=1 ["rgg_n_2_19_s0"]=1 ["soc-delicious"]=1 ["ca-coauthors-dblp"]=1
    ["soc-digg"]=1 ["eu-2005"]=1 ["web-Google"]=1 
    ["rgg_n_2_20_s0"]=1 ["soc-pokec"]=1 
    ["rgg_n_2_21_s0"]=1 ["rgg_n_2_22_s0"]=1
    ["rgg_n_2_23_s0"]=1 ["hugetrace-00010"]=1 ["hugetrace-00020"]=1 ["rgg_n_2_24_s0"]=1
)

while read -r line; do
    [[ "$line" == "*********"* ]] && continue
    read -r graph_file _ _ _ _ density _ <<< "$line"
    if [[ "$density" == "-nan" ]]; then
        continue
    fi

    clean_path="${graph_file#../../}"
    clean_path="${clean_path//\/\//\/}"
    
    graph_name=$(basename "$clean_path")

    if [[ ! -v allowed_graphs[$graph_name] ]]; then
        continue
    fi
    
    truncated_density=$(echo "$density" | awk '{
        split($0,a,".");
        if (length(a[2])>=2) printf "%s.%s", a[1], substr(a[2],1,2);
        else printf "%s", $0
    }')

    for seed in {1..10}; do
    (
        result_line=$("$NUQCLQ_PATH" "$clean_path" "$truncated_density" 2 0.0001 "$seed" 2>/dev/null)

        if [[ -n "$result_line" ]]; then
            echo "$result_line" >> "edqc-result-seed-$seed.out"
        fi
    ) &
    done
    wait
done < "$RESULTS_FILE"