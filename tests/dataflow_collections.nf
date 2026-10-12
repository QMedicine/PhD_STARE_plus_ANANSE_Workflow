nextflow.enable.dsl=2
workflow {
    networks = channel.of(tuple('historic','Old','old.tsv'),tuple('historic','Young','young.tsv')).collect(flat:false).map { nets -> [nets] }
    des = channel.of(tuple('Old','Young','forward.tsv'),tuple('Young','Old','reverse.tsv')).collect(flat:false).map { ds -> [ds] }
    channel.of(tuple('Old','Young'),tuple('Young','Old')).combine(networks).combine(des).map { source,target,nets,ds ->
        assert nets.size() == 2
        assert ds.size() == 2
        tuple(source,target,nets.find { n -> n[1]==source }[2],nets.find { n -> n[1]==target }[2],ds.find { d -> d[0]==source && d[1]==target }[2])
    }.view()
}
