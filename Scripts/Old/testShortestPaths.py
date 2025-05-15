from gothic import *


if __name__ == "__main__":
    # for each case
    for case in ["MEC1c", "MEC2a", "MEC2b", "MEC2c"]:
        for flevel in [0.05, 0.1, 0.2]:

            # initialize vars
            fpstring = "Top" + str(int(flevel * 100)) + "p"

            prefix = case + FixedBins500kb
            norm_graph = prefix + "_ICEnormed.gt"
            tf_graph = prefix + "_ICEnormedEDGEMOD.gt"
            filt_norm_graph = prefix + "_ICEnormed" + fpstring + ".gt"
            filt_tf_graph = prefix + "_ICEnormedEDGEMOD" + fpstring + ".gt"
            norm_distmat = prefix + "_ICEnormed" + fpstring + "_APSP.tsv"
            tf_distmat = prefix + "_ICEnormedEDGEMOD" + fpstring + "_APSP.tsv"
            norm_edge_plot = prefix + "_ICEnormed" + fpstring + "_edgeplot.png"
            norm_apsp_plot = prefix + "_ICEnormed" + fpstring + "_APSPplot.png"
            tf_edge_plot = prefix + "_ICEnormedEDGEMOD" + fpstring + "_edgeplot.png"
            tf_apsp_plot = prefix + "_ICEnormedEDGEMOD" + fpstring + "_APSPplot.png"
            norm_report = prefix + "_ICEnormed" + fpstring + "REPORT.txt"
            tf_report = prefix + "_ICEnormedEDGEMOD" + fpstring + "REPORT.txt"

            # transform weights
            modify_weights(norm_graph, outfile=tf_graph)

            # filter
            create_top_edges_graph(norm_graph, threshold=flevel, outfile=filt_norm_graph, weights_are_distances=False)
            create_top_edges_graph(tf_graph, threshold=flevel, outfile=filt_tf_graph, weights_are_distances=True)

            # get APSP
            annotate_apsp(filt_norm_graph, filt_norm_graph)
            distmatrix_from_annotations(filt_norm_graph, outfile=norm_distmat)
            annotate_apsp(filt_tf_graph, filt_tf_graph)
            distmatrix_from_annotations(filt_tf_graph, outfile=tf_distmat)

            # plot
            plot_edge_weight_distribution(filt_norm_graph, outfile=norm_edge_plot)
            plot_shortest_path_pairs_distribution(filt_norm_graph, outfile=norm_edge_plot)
            plot_edge_weight_distribution(filt_tf_graph, outfile=norm_tf_plot)
            plot_shortest_path_pairs_distribution(filt_tf_graph, outfile=norm_tf_plot)

            # generate graph reports
            generate_graph_report(filt_norm_graph)
            generate_graph_report(filt_tf_graph)
