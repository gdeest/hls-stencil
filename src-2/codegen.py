#!/usr/bin/env python


from islpy import *

def separator():
    print("")
    # print("-------------------------------------------------------------------------------------------")

## Builds the following map: {x -> x' | x' = x . size + offset }

def build_affine_expansion(space, size, offset):
    n = space.dim(dim_type.set)

    map_space = Space.map_from_domain_and_range(space, space)
    lspace = LocalSpace.from_space(map_space)
    expansion = Map.universe(map_space)

    for i, (s_i, o_i) in enumerate(zip(size, offset)):
        # x_i' - s_i * x_i - o_i = 0
        constraint = Constraint.equality_alloc(lspace)

        constraint = constraint.set_coefficient_val(dim_type.in_, i, s_i)
        constraint = constraint.set_coefficient_val(dim_type.out, i, -1)
        constraint = constraint.set_constant_val(o_i)

        expansion = expansion.add_constraint(constraint)

    return expansion

## Builds the following maps:
#
#  - dense_to_origins: Maps points from a dense space to a regular lattice of tile origins such that offset is the image of zero.
#    More formally: { dense[x] -> origins[x'] | x' = x . size + offset }
#
#  - origins_to_points { origins[x] -> D[x'] | (exists k. x = k . x + offset) and (x <= x' < x + size) }
#    Maps points tile origins to points within the tile.
#
#  Optionally, names can be specified for the dense space and the space of tile origins. Otherwise, the original name is kept.

def tile_space(space, size, offset, origins_lattice_name=None, dense_lattice_name=None):
    if origins_lattice_name is None:
        origins_lattice_name = space.get_tuple_name(dim_type.set)
    if dense_lattice_name is None:
        dense_lattice_name = origins_lattice_name

    n = space.dim(dim_type.set)

    dense_space = space.set_tuple_name(dim_type.set, dense_lattice_name)
    dense_to_origins = build_affine_expansion(dense_space, size, offset)
    dense_to_origins = dense_to_origins.set_tuple_name(dim_type.out, origins_lattice_name)

    origins_space = dense_to_origins.get_space().range()
    origins_to_points_space = Space.map_from_domain_and_range(origins_space, space)

    origins_to_points = Map.universe(origins_to_points_space)
    otp_lspace = LocalSpace.from_space(origins_to_points_space)

    for i, (s_i, o_i) in enumerate(zip(size, offset)):
        # x'i - xi >= 0
        constraint = Constraint.inequality_alloc(otp_lspace)
        constraint = constraint.set_coefficient_val(dim_type.out, i, 1)
        constraint = constraint.set_coefficient_val(dim_type.in_, i, -1)
        origins_to_points = origins_to_points.add_constraint(constraint)

        # x_i' < x_i + s  ->  x_i - x_i' + (s-1) >= 0
        constraint = Constraint.inequality_alloc(otp_lspace)
        constraint = constraint.set_coefficient_val(dim_type.out, i, -1)
        constraint = constraint.set_coefficient_val(dim_type.in_, i, 1)
        constraint = constraint.set_constant_val(s_i-1)
        origins_to_points = origins_to_points.add_constraint(constraint)

    origins_to_points = origins_to_points.intersect_domain(dense_to_origins.range())

    #return (dense_to_origins, origins_to_points)
    return (dense_to_origins, origins_to_points)

def remove_self(map):
    # Yet again this stupid bug.
    map = Map("%s" % map)
    return map.subtract(Map.identity(map.get_space()))

def point_to_tuple(point, dimtype = dim_type.set):
    tuple = ()
    for idx in range(point.get_space().dim(dimtype)):
        tuple = tuple + (point.get_coordinate_val(dimtype, idx).get_num_si(),)
    return tuple

def compute_dep_vectors(deps):
    # Stupid isl bug workaround
    deps = Map("%s" % deps)

    # Remove any parameter
    n = deps.domain().dim(dim_type.param)
    deps = deps.remove_dims(dim_type.param, 0, n)

    # Create origin in parameter-less space
    zero = Point.zero(deps.domain().get_space())
    zero = Set.from_point(zero)

    # Create vector list
    vecs = []
    deps.intersect_domain(zero).range().foreach_point(lambda p: vecs.append(point_to_tuple(p)))

    return vecs

def compute_inner_outer_tiles(D, R, deps, dense_to_points):
    assert(R.is_subset(D))

    # Compute the /interior/ region of D, points which:
    # - are in R
    # - do not depend from points outside R

    not_in_I = deps.intersect_domain(R).subtract_range(R).domain()
    I = R.subtract(not_in_I)

    # Set of all non-empty tiles.
    all_tiles   = dense_to_points.intersect_range(D).domain()
    # Set of tiles all tiles (be they empty or not) having points /outside/ of R.
    non_regular = dense_to_points.subtract_range(R).domain()
    # Regular tiles: tiles entirely contained in the regular region.
    regular     = all_tiles.subtract(non_regular)
    # Tiles which have dependencies outside of D.
    out_deps = dense_to_points.apply_range(deps).subtract_range(D).domain()
    # Inner tiles are regular tiles which have no dependencies outside of D.
    inner       = regular.subtract(out_deps)
    outer       = all_tiles.subtract(inner)

    return(all_tiles, inner, outer)

def compute(D, R, deps, skew, tile_skew, tile_size, tile_offset):
    assert(skew.get_tuple_name(dim_type.in_) == "D")
    assert(skew.get_tuple_name(dim_type.out) == "D_s")

    D_s = skew.intersect_domain(D).range()
    R_s = skew.intersect_domain(R).range()

    skewed_deps = skew.reverse().apply_range(deps).apply_range(skew)
    skewed_2 = deps.apply_domain(skew).apply_range(skew)

    (dto, otp) = tile_space(skew.get_space().range(), tile_size, tile_offset, "O_T", "T")
    dense_to_points = dto.apply_range(otp)
    skewed_to_points = tile_skew.reverse().apply_range(dense_to_points)
    skewed_tiles_to_orig = skewed_to_points.apply_range(skew.reverse())

    (all_tiles, inner, outer) = compute_inner_outer_tiles(D_s, R_s, skewed_deps, dense_to_points)

    n_dim = all_tiles.n_dim()

    wavefronts = (tile_skew
                  .intersect_domain(all_tiles)
                  .range()
                  .project_out(dim_type.set, 1, n_dim-1)
                  .set_tuple_name("WF"))
    # print("Wavefronts: ", wavefronts)

    wf_space = wavefronts.get_space()
    wavefront_schedule_space = Space.map_from_domain_and_range(wf_space, wf_space)
    wavefront_schedule = Map.identity(wavefront_schedule_space).intersect_domain(wavefronts)
    # print("Wavefront schedule: ", wavefront_schedule)

    sched = tile_skew.move_dims(dim_type.param, 0, dim_type.out, 0, 1).set_dim_name(dim_type.param, 0, "d")
    inner_sched = sched.intersect_domain(inner)
    outer_sched = sched.intersect_domain(outer)

    # print("Per-wavefront inner schedule: ", outer_sched)
    # print("Per-wavefront outer schedule: ", inner_sched)

    build = AstBuild.from_context(Set("{:}"))
    ast = build.ast_from_schedule(wavefront_schedule)

    idT = sched.domain().identity()
    inner_sched = inner_sched.apply_domain(idT.set_tuple_name(dim_type.out, "INNER_TILE"))
    outer_sched = outer_sched.apply_domain(idT.set_tuple_name(dim_type.out, "OUTER_TILE"))

    inner_ast = build.ast_from_schedule(inner_sched)
    outer_ast = build.ast_from_schedule(outer_sched)

    def callback(ast, p, options, node, thread_name, before_wf=[], after_wf=[]):
        p = p.start_line()
        p = p.print_str("{\n")
        p = p.set_prefix("\t")

        p = p.start_line()
        d = node.user_get_expr().get_op_arg(1)
        p = p.print_str("int d = ")
        p = p.print_ast_expr(d)
        p = p.print_str(";\n")

        p = p.start_line()
        p = p.print_str("sem_wait(mutex_cout);\n")
        p = p.start_line()
        p = p.print_str("cout << \"%s blocked on barrier: \" << (d&1) << endl;\n" % thread_name)
        p = p.print_str("sem_post(mutex_cout);\n")
        p = p.start_line()
        p = p.print_str("pthread_barrier_wait(&barriers[d & 1]);\n")
        p = p.start_line()
        p = p.print_str("sem_wait(mutex_cout);\n")
        p = p.start_line()
        p = p.print_str("cout << \"%s unblocked.\" << endl;\n" % thread_name)
        p = p.start_line()
        p = p.print_str("sem_post(mutex_cout);\n")

        p = p.start_line()
        p = p.print_str("// Before\n")
        for line in before_wf:
            p = p.start_line()
            p = p.print_str(line + "\n")
        p = p.start_line()
        p = p.print_str("\n")

        pp = Printer.to_str(DEFAULT_CONTEXT)
        print_options = AstPrintOptions.alloc(DEFAULT_CONTEXT)
        pp = ast.print_(pp, print_options)
        for line in pp.get_str().split('\n'):
            p = p.start_line()
            p = p.print_str(line + "\n")

        p = p.start_line()
        p = p.print_str("// After WF\n")
        for line in after_wf:
            p = p.start_line()
            p = p.print_str(line + "\n")
        p = p.start_line()
        p = p.print_str("\n")

        p = p.set_prefix("")
        p = p.start_line()
        p = p.print_str("}\n")

        return p

    before_wf_inner = [ "int ntiles = 0;",
                        "std::vector<int> coords;"]
    size_params = ", ".join(["T", "N"])
    after_wf_inner  = [ "if (ntiles > 0) {",
                        "  int *coords_ = (int *) ALLOC(%d * ntiles * sizeof(int));" % n_dim,
                        "  for (int i=0; i<%d*ntiles; i++)" % n_dim,
                        "    coords_[i] = coords[i];",
                        "  inner_tiles_wrapper(arr, %s, ntiles, coords_, mutex_cout, cnt);" % size_params,
                        "  FREE(coords_);",
                        "}"]

    inner_callback = (lambda p, options, node: callback(inner_ast, p, options, node, "HW", before_wf_inner, after_wf_inner))
    outer_callback = (lambda p, options, node: callback(outer_ast, p, options, node, "SW"))


    def wf_callback(p, options, node):
        lines = [
            "pthread_barrier_wait(&barriers[d&1]);",
            "sem_wait(mutex_cout);",
            "cout << \"Wavefront: \" << d << endl;",
            "sem_post(mutex_cout);",
            "pthread_barrier_destroy(&barriers[d&1]);",
            "pthread_barrier_init(&barriers[d&1], 0, 3);"
        ]

        p = p.start_line()
        p = p.print_str("{\n")
        p = p.set_prefix("  ")

        p = p.start_line()
        d = node.user_get_expr().get_op_arg(1)
        p = p.print_str("int d = ")
        p = p.print_ast_expr(d)
        p = p.print_str(";\n")

        for l in lines:
            p = p.start_line()
            p = p.print_str("%s\n" % l)

        p = p.set_prefix("")
        p = p.start_line()
        p = p.print_str("}\n")

        return p


    def to_str(cb):
        print_options = AstPrintOptions.alloc(DEFAULT_CONTEXT)
        print_options, cb_handle = print_options.set_print_user(cb)
        p = Printer.to_str(DEFAULT_CONTEXT)
        p = p.set_output_format(format.C)

        p = ast.print_(p, print_options)

        return p.get_str()

    return (to_str(inner_callback), to_str(outer_callback), to_str(wf_callback))


    # wavefronts = tile_skew.intersect_domain(all_tiles).range().project_out(dim_type.set, 1, n_dim-1)
    # all_tiles_skewed = tile_skew.intersect_domain(all_tiles).range()
    # inner_skewed = tile_skew.intersect_domain(inner).range()
    # outer_skewed = tile_skew.intersect_domain(outer).range()

    # print("All tiles:\t", all_tiles)
    # print("Inner_tiles:\t", inner)
    # print("Outer_tiles:\t", outer)
    # separator()

    # tile_deps = remove_self(dense_to_points.apply_range(skewed_deps).apply_range(dense_to_points.reverse()))
    # skewed_tile_deps = tile_skew.reverse().apply_range(tile_deps).apply_range(tile_skew)

    # print("Inter-tile dependencies (after tile skew): ")
    # dep_vectors = compute_dep_vectors(skewed_tile_deps)
    # for vec in dep_vectors:
    #     print(vec)
    # separator()

    # n_dim = all_tiles.n_dim()

    # #wavefronts = all_tiles.project_out(dim_type.set, 1, n_dim-1)
    # range_space  = Space.set_alloc(DEFAULT_CONTEXT, 0,0)
    # domain_space = Space.set_alloc(DEFAULT_CONTEXT, 0, n_dim)
    # map_space = Space.map_from_domain_and_range(domain_space, range_space)
    # universe = Map.universe(map_space)
    # universe = Map("%s" % universe)

    # def wavefronts_to_tiles(tile_set):
    #     sets = []
    #     tile_set.foreach_set(sets.append)
    #     sets = list(map(lambda s: s.set_tuple_name(""), sets))
    #     tile_set = UnionSet("{}")
    #     for set in sets:
    #         tile_set = tile_set.union(set)

    #     tile_set = Set("%s" % tile_set)
    #     m = universe.intersect_domain(tile_set)
    #     m = m.move_dims(dim_type.out, 0, dim_type.in_, 1, n_dim-1)

    #     m = m.insert_dims(dim_type.out, 0, 1)


    #     return m

    # #print(wavefronts_to_tiles(all_tiles_skewed))
    # print("Inner wavefronts: ", wavefronts_to_tiles(inner_skewed))
    # print("Outer wavefronts: ", wavefronts_to_tiles(outer_skewed))

def build_space(ndims):
    space = Space.set_alloc(DEFAULT_CONTEXT, 2, ndims+1)
    space = space.set_tuple_name(dim_type.set, "D")

    space = space.set_dim_name(dim_type.param, 0, "T")
    space = space.set_dim_name(dim_type.param, 1, "N")
    # for i in range(1, ndims+1):
    #     space = space.set_dim_name(dim_type.param, i, "N%d" % (i-1))

    return space

def build_domain(ndims):
    if ndims < 1:
        raise Exception("Error: number of dimensions must be at least 1.")

    space = build_space(ndims)

    domain = BasicSet.universe(space)
    local_space = LocalSpace.from_space(space)

    for i in range(0, ndims+1):
        constraint = Constraint.alloc_inequality(local_space)

        constraint = constraint.set_coefficient_val(dim_type.set, i, 1)
        constraint = constraint.set_constant_val(-1)
        domain = domain.add_constraint(constraint)

        param_dim = 1
        if i == 0:
            param_dim = 0

        constraint = Constraint.alloc_inequality(local_space)
        constraint = constraint.set_coefficient_val(dim_type.set, i, -1)
        constraint = constraint.set_coefficient_val(dim_type.param, param_dim, 1)
        if i>0:
            constraint = constraint.set_constant_val(-2)
        else:
            constraint = constraint.set_constant_val(-1)
        domain = domain.add_constraint(constraint)

    return domain

def build_skew(ndims, mat):
    space = Space.alloc(DEFAULT_CONTEXT, 0, ndims, ndims)
    local_space = LocalSpace.from_space(space)

    skew = Map.universe(space)
    for i in range(ndims):
        constraint = Constraint.alloc_equality(local_space)

        for j in range(ndims):
            constraint = constraint.set_coefficient_val(dim_type.in_, j, mat[i][j])

        constraint = constraint.set_coefficient_val(dim_type.out, i, -1)
        skew = skew.add_constraint(constraint)

    return skew

def build_dependence_relation(ndims, dep_vectors):
    space = build_space(ndims)

    space = Space.map_from_domain_and_range(space, space)
    local_space = LocalSpace.from_space(space)

    umap = UnionMap.empty(space)
    for v in dep_vectors:
        map = Map.universe(space)

        for i in range(ndims+1):
            constraint = Constraint.alloc_equality(local_space)
            constraint = constraint.set_coefficient_val(dim_type.in_, i, -1)
            constraint = constraint.set_coefficient_val(dim_type.out, i,  1)
            constraint = constraint.set_constant_val(-v[i])
            map = map.add_constraint(constraint)

        umap = umap.add_map(map)

    return umap

def parse_matrix(ndims, str):
    rows = str.split(";")
    if len(rows) != ndims:
        raise Exception("Error: invalid number of rows.")

    mat = []
    for r in rows:
        elts = list(map(lambda s: int(s.strip()), r.split(",")))
        if len(elts) != ndims:
            raise Exception("Error: row has an invalid number of elements.")
        mat.append(elts)

    return mat

def parse_tuple(str):
    return tuple(map(lambda s: int(s.strip()), str.split(",")))

if __name__ == '__main__':
    import sys
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-o", "--output-file", action="store", type="string", dest="output_file")
    parser.add_option("-d", "--data-dims", action="store", type="int", dest="ndims")
    # parser.add_option("-R", "--regular-predicate", action="store", type="string", dest="regular")
    parser.add_option("-S", "--tile-size", action="store", type="string", dest="tile_size")
    parser.add_option("--skew", action="store", type="string", dest="skew")
    parser.add_option("--tile-skew", action="store", type="string", dest="tile_skew")
    parser.add_option("--dep", action="append", type="string", dest="deps", default=[])

    (options, args) = parser.parse_args(sys.argv)

    # separator()
    if options.ndims is None:
        parser.error("Number of data dimensions (-d, --data-dims) is required.")
    else:
        ndims = options.ndims

    D = build_domain(ndims)
    print("Domain: ", D)

    R = D

    # if options.regular is not None:
    #     R = Set("[%s] -> {[%s]: %s}" % (",".join(build_params(ndims), build_indices(ndims)), options.regular))
    # else:
    #     R = D
    # print("Regular region: ", R)

    if options.tile_size is None:
        parser.error("Tile size is require (-S / --tile-size)")
    else:
        tile_size = parse_tuple(options.tile_size)

    if options.skew is None:
        parser.error("Skewing required (--skew)")
    else:
        skew = build_skew(ndims+1, parse_matrix(ndims+1, options.skew))
        skew = skew.set_tuple_name(dim_type.in_, "D")
        skew = skew.set_tuple_name(dim_type.out, "D_s")

    if options.tile_skew is None:
        parser.error("Tile skewing required (--tile-skew)")
    else:
        tile_skew = build_skew(ndims+1, parse_matrix(ndims+1, options.tile_skew))
        tile_skew = tile_skew.set_tuple_name(dim_type.in_, "T")
        tile_skew = tile_skew.set_tuple_name(dim_type.out, "T_s")

    # print("Skew: ", skew)

    dep_vectors = list(map(parse_tuple, options.deps))
    deps = build_dependence_relation(ndims, dep_vectors)
    # print("Dependence relation: ", deps)
    # separator()

    tile_offset = ()
    for i in range(ndims+1):
        tile_offset = tile_offset + (0,)

    inner_code, outer_code, wf_code = compute(D, R, deps, skew, tile_skew, tile_size, tile_offset)
    template_params = (",".join(["%d"] * (ndims+1))) % tile_size
    # domain_args = ",".join(["int T"] + ["int N%d" % i for i in range(0, ndims)])
    domain_args = ",".join(["int T", "int N"])
    macro_params = ", ".join(["t"] + ["x%d" % i for i in range(0, ndims)])
    wrapper_params = ", ".join(["T", "N"]) + ", "
    # wrapper_params = ", ".join(["T"] + ["N%d" % i for i in range(0, ndims)]) + ", "
    wrapper_params = wrapper_params + ", ".join(["(t)"] + ["(x%d)" % i for i in range(0, ndims)])

    coords_push = []
    for coord in (["t"] + ["x%d" % i for i in range(0, ndims)]):
        coords_push.append("    coords.push_back(%s); \\" % coord)
    coords_push = '\n'.join(coords_push)

    inner_code = '\n'.join([
        # "#define INNER_TILE(%s) inner_tile_wrapper(arr, %s, mutex_cout, cnt)" % (macro_params, wrapper_params),
        "#define INNER_TILE(%s) \\" %macro_params,
        "  {\\",
        coords_push,
        "    ntiles++; \\",
        "  }",
        "template<> void scan_hw_tiles<%s>(%s, data_t *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt) {" % (template_params, domain_args),
        inner_code,
        "}\n",
        "#undef INNER_TILE"
    ])

    outer_code = '\n'.join([
        "#define OUTER_TILE(%s) outer_tile_wrapper(arr, %s, mutex_cout, cnt)" % (macro_params, wrapper_params),
        "template<> void scan_sw_tiles<%s>(%s, data_t *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt) {" % (template_params, domain_args),
        outer_code,
        "}\n",
        "#undef OUTER_TILE"
    ])

    wf_code = '\n'.join([
        "template<> void scan_wavefronts<%s>(%s, pthread_barrier_t barriers[2], sem_t *mutex_cout) {" % (template_params, domain_args),
        wf_code,
        "}\n"
    ])

    print("#undef max")
    print("#undef min")
    print("#include <vector>")
    print(inner_code)
    print(outer_code)
    print(wf_code)
