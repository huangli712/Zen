using Documenter

makedocs(
    sitename="Flink: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/Flink",
        size_threshold = 409600, # 400kb
        assets = ["assets/flink.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = Module[],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => [
            "About Flink" => "intro.md",
            "Installation" => "install.md",
            "Basic Usage" => "usage.md",
        ],
        "Modules" => Any[
            "Constants" => "guide/m_constants.md",
            "Colorful Terminal Output" => "guide/m_face.md",
            "Linked List" => "guide/m_linkedlist.md",
            "Message Passing Interface" => "guide/m_mpi.md",
            "Configuration Parser" => "guide/m_parser.md",
            "Sparse Matrix" => "guide/m_sparse.md",
            "Pseudorandom Number Generator" => "guide/m_spring.md",
            "Stack" => "guide/m_stack.md",
            "Analytical Tetrahedron Algorithm" => "guide/m_tetra.md",
        ],
        "Subroutines" => Any[
            "Error and Exception" => "guide/s_error.md",
            "Fourier Transformation" => "guide/s_fourier.md",
            "Special Functions" => "guide/s_function.md",
            "Integration" => "guide/s_integrator.md",
            "Matrix" => "guide/s_matrix.md",
            "Spline Interpolation" => "guide/s_spline.md",
            "Utility" => "guide/s_util.md",
            "Vector" => "guide/s_vector.md",
        ],
        "Appendix" => Any[
            "make.inc" => "appendix/make.md",
        ],
    ],
)
