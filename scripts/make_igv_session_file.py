import sys, glob, os, re, argparse
import matplotlib.colors as mcolors
import pandas as pd
    
def make_igv_session_file(bigwigs_fnames, to_group_function, out_file_name="igv_session.from_script.xml"):
    fnames = sorted(list(bigwigs_fnames))

    groups = sorted(list(set([to_group_function(x) for x in fnames])))

    skip_colors = [f"xkcd:{cn}" for cn in ["bland", "booger", "electric lime", "yellowish tan", "really light blue", "cloudy blue", "egg shell", "bright seq green"]]
    
    color_table = {k: ','.join([str(int(255*x)) for x in mcolors.to_rgb(c)]) \
                   for k,c in mcolors.XKCD_COLORS.items() if k not in skip_colors}
        
    colors_list = list(color_table.values())
    n_colors = len(list(color_table.values()))
    n_groups = len(groups)


    def to_color(fname):
        return colors_list[groups.index(to_group(fname)) % n_colors]
    
    header = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>

    <Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" locus="chr5:6612192-6649998" nextAutoscaleGroup="13" version="8">

        <Resources>

    """
    resource_lines = '\n'.join([f"""
            <Resource path="{os.path.realpath(fname)}" type="bigwig"/>""" for fname in fnames])

    header2 = """

        </Resources>

        <Panel height="2480" name="DataPanel" width="1625">

        """

    track_lines = '\n'.join([f"""    <Track attributeKey="{os.path.basename(fname)}" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" color="{to_color(fname)}" fontSize="10" id="{os.path.realpath(fname)}" name="{os.path.basename(fname)}" renderer="BAR_CHART" visible="true" windowFunction="mean">

                <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="51.3913" minimum="0.0" type="LINEAR"/>

            </Track>""" for fname in fnames])

    footer = """
        </Panel>

        <Panel height="299" name="FeaturePanel" width="1625">

            <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="POSITIVE" shouldShowTranslation="false" visible="true"/>

            <Track attributeKey="Gene" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;770.0;255,255,255;0,0,178" displayMode="EXPANDED" fontSize="10" groupByStrand="false" height="35" id="hg38_genes" name="Gene" visible="true"/>

        </Panel>

        <PanelLayout dividerFractions="0.6633444075304541"/>

        <Regions>

            <Region chromosome="chr12" end="53431282" start="53418824"/>

        </Regions>

        <HiddenAttributes>

            <Attribute name="DATA FILE"/>

            <Attribute name="DATA TYPE"/>

            <Attribute name="NAME"/>

        </HiddenAttributes>

    </Session>
    """

    outstr = header + resource_lines + header2 + track_lines + footer

    with open(out_file_name, 'w') as f:
        f.write(outstr)
        
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--signal_dir", help="Bigwig folder.")
    parser.add_argument("--samplesheet", help="Samplesheet.")
    parser.add_argument("--outfilename", help="Output igv session xml filename.")
    
    args = parser.parse_args()

    bigwigs_fnames = glob.glob(f"{args.signal_dir}/*scaled.bigwig")
    
    # Use the samplesheet to associate samples with conditions for coloring.
    ss = pd.read_csv(args.samplesheet).set_index('sample')
    
    def to_group(bwfname):
        return ss.loc[re.search('/([^/]+).scaled.bigwig', bwfname).group(1), 'condition']

    make_igv_session_file(bigwigs_fnames, to_group, out_file_name=args.outfilename)



