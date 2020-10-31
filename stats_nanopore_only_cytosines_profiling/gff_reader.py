
def extract_region_by_feature(gff3_eles, gfeature):
    outstrs = []
    for gff_ele in gff3_eles:
        if gff_ele.get_feature() == gfeature:
            outstrs.append(gff_ele)
    # print("extract {} regions by feature ({})".format(len(outstrs), gfeature))
    return outstrs


def extract_region_by_features(gff3_eles, gfeatures):
    outstrs = []
    gfeatures = set(gfeatures)
    for gff_ele in gff3_eles:
        if gff_ele.get_feature() in gfeatures:
            outstrs.append(gff_ele)
    print("extract {} regions by feature ({})".format(len(outstrs), gfeatures))
    return outstrs


def extract_region_by_note(gff3_eles, gnote):
    outstrs = []
    for gff_ele in gff3_eles:
        if "Note" in gff_ele.get_attr_keys() and gff_ele.get_attrs()["Note"] == gnote:
            outstrs.append(gff_ele)
    print("extract {} regions by Note ({})".format(len(outstrs), gnote))
    return outstrs


def extract_region_by_notes(gff3_eles, gnotes):
    gnotes = set(gnotes)
    outstrs = []
    for gff_ele in gff3_eles:
        if "Note" in gff_ele.get_attr_keys() and gff_ele.get_attrs()["Note"] in gnotes:
            outstrs.append(gff_ele)
    print("extract {} regions by Notes ({})".format(len(outstrs), gnotes))
    return outstrs


def extract_region_by_biotype(gff3_eles, gtype):
    outstrs = []
    for gff_ele in gff3_eles:
        if "gene_biotype" in gff_ele.get_attr_keys() and gff_ele.get_attrs()["gene_biotype"] == gtype:
            outstrs.append(gff_ele)
    print("extract {} regions by gene_biotype ({})".format(len(outstrs), gtype))
    return outstrs


def extract_region_by_biotypes(gff3_eles, gtypes):
    gtypes = set(gtypes)
    outstrs = []
    for gff_ele in gff3_eles:
        if "gene_biotype" in gff_ele.get_attr_keys() and gff_ele.get_attrs()["gene_biotype"] in gtypes:
            outstrs.append(gff_ele)
    print("extract {} regions by gene_biotypes ({})".format(len(outstrs), gtypes))
    return outstrs


def extract_region_by_attri(gff3_eles, attri_name, attri_val):
    outstrs = []
    for gff_ele in gff3_eles:
        if attri_name in gff_ele.get_attr_keys() and gff_ele.get_attrs()[attri_name] == attri_val:
            outstrs.append(gff_ele)
    # print("extract {} regions by {} ({})".format(len(outstrs), attri_name, attri_val))
    return outstrs


def extract_region_by_attri_conditionFeature(gff3_eles, attri_name, attri_val, feature):
    outstrs = []
    for gff_ele in gff3_eles:
        if gff_ele.get_feature() == feature:
            if attri_name in gff_ele.get_attr_keys() and gff_ele.get_attrs()[attri_name] == attri_val:
                outstrs.append(gff_ele)
    # print("extract {} regions by {} ({})".format(len(outstrs), attri_name, attri_val))
    return outstrs


def extract_region_by_notattri_conditionFeature(gff3_eles, attri_name, attri_val, feature):
    outstrs = []
    for gff_ele in gff3_eles:
        if gff_ele.get_feature() == feature:
            if attri_name in gff_ele.get_attr_keys() and gff_ele.get_attrs()[attri_name] != attri_val:
                outstrs.append(gff_ele)
    # print("extract {} regions by {} ({})".format(len(outstrs), attri_name, attri_val))
    return outstrs


# ===================================
def get_kinds_of_a_attri(gff3_eles, attri_name):
    attri_kinds = set()
    for gff_ele in gff3_eles:
        if attri_name in gff_ele.get_attr_keys():
            attri_kinds.add(gff_ele.get_attrs()[attri_name])
    return attri_kinds


def get_kinds_of_a_attri_conditionFeature(gff3_eles, attri_name, feature):
    attri_kinds = set()
    for gff_ele in gff3_eles:
        if gff_ele.get_feature() == feature:
            if attri_name in gff_ele.get_attr_keys():
                attri_kinds.add(gff_ele.get_attrs()[attri_name])
    return attri_kinds


def get_featurekinds_of_a_attrikey(gff3_eles, attri_key):
    feature_kinds = set()
    for gff_ele in gff3_eles:
        if attri_key in gff_ele.get_attr_keys():
            feature_kinds.add(gff_ele.get_feature())
    return feature_kinds


def get_featurekinds_conditionIDstarts(gff3_eles, IDstartswith):
    feature_kinds = set()
    for gff_ele in gff3_eles:
        if gff_ele.get_id().startswith(IDstartswith):
            feature_kinds.add(gff_ele.get_feature())
    return feature_kinds


class GFF3Element:
    def __init__(self, fields):
        # "chromosome", "source", "feature", "start", "end",
        #              "score", "strand", "phase", "attributes"
        self._chromosome = fields[0]
        self._source = fields[1]
        self._feature = fields[2]
        self._start = int(fields[3]) - 1  # turn to 0-based
        self._end = int(fields[4])
        self._score = fields[5]
        self._strand = fields[6]
        self._phase = fields[7]
        self._attributes = fields[8]

        self._set_gene_attrs()

    def _set_gene_attrs(self):
        self._attrs = dict()
        for attr_kv in self._attributes.strip().split(";"):
            if attr_kv != "":
                attr = attr_kv.strip().split("=")
                self._attrs[attr[0]] = attr[1]
        self._attr_keys = set(self._attrs.keys())

    def get_chromosome(self):
        return self._chromosome

    def get_source(self):
        return self._source

    def get_feature(self):
        return self._feature

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    def get_score(self):
        return self._score

    def get_strand(self):
        return self._strand

    def get_phase(self):
        return self._phase

    def get_attributes(self):
        return self._attributes

    def get_attrs(self):
        return self._attrs

    def get_attr_keys(self):
        return self._attr_keys

    def get_id(self):
        if "ID" in self._attr_keys:
            return self._attrs["ID"]
        elif "Name" in self._attr_keys:
            return self._attrs["Name"]
        return "-"

    def print_str(self):
        # id, chrom, start(0-based), end(0-based), strand, feature, attributes
        return "\t".join([self.get_id(), self.get_chromosome(), str(self.get_start()),
                          str(self.get_end()), self.get_strand(), self.get_feature(),
                          self.get_attributes()])


class GFF3:
    def __init__(self, filepath):
        self.eles = []
        self._features = set()
        with open(filepath, "r") as rf:
            for line in rf:
                if not line.startswith("#"):
                    words = line.strip().split("\t")
                    gffele = GFF3Element(words)
                    self.eles.append(gffele)
                    self._features.add(gffele.get_feature())

    def get_eles(self):
        return self.eles

    def get_features(self):
        return self._features
