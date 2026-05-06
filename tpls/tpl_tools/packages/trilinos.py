from tpl_tools.packages import packages
from tpl_tools import utils
import os
import sys

# https://patch-diff.githubusercontent.com/raw/trilinos/Trilinos/pull/15159.patch
sacado_patch = """
diff --git a/packages/sacado/src/Kokkos_LayoutContiguous.hpp b/packages/sacado/src/Kokkos_LayoutContiguous.hpp
index 316b7796749..5c905ba17e9 100644
--- a/packages/sacado/src/Kokkos_LayoutContiguous.hpp
+++ b/packages/sacado/src/Kokkos_LayoutContiguous.hpp
@@ -18,6 +18,7 @@
 #endif
 #include "Kokkos_Core_fwd.hpp"
 #include "Kokkos_Layout.hpp"
+#include "Kokkos_DynRankView.hpp"
 #ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
 #undef KOKKOS_IMPL_PUBLIC_INCLUDE
 #undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
@@ -83,36 +84,116 @@ struct inner_layout< LayoutContiguous<Layout, Stride> > {
   typedef Layout type;
 };
 
-} // namespace Kokkos
-
-// FIXME This is evil and needs refactoring urgently.
-// Make LayoutContiguous<Layout> equivalent to Layout
-namespace std {
-
-  template <class Layout, unsigned Stride>
-  struct is_same< Kokkos::LayoutContiguous<Layout,Stride>, Layout> {
-    static const bool value = true;
-  };
+namespace Impl {
 
-  template <class Layout, unsigned Stride>
-#if defined(KOKKOS_COMPILER_INTEL)
-  inline constexpr bool is_same_v< Kokkos::LayoutContiguous<Layout,Stride>, Layout> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
+  // Specialize DynRankDimTraits for LayoutContiguous.
+  // Weirdly, we need a full specialization on bool for the non-legacy View impl case
+  // because of this code in Kokkos_DynRankView.hpp (around line 429):
+  // #ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
+  //   using drdtraits = Impl::DynRankDimTraits<typename view_type::specialize>;
+  // #else
+  //   using drdtraits = Impl::DynRankDimTraits<
+  //       std::conditional_t<view_type::traits::impl_is_customized, bool, void>>;
+  // #endif
+  template <>
+#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
+  struct DynRankDimTraits<ViewSpecializeSacadoFadContiguous> {
 #else
-  static constexpr bool is_same_v< Kokkos::LayoutContiguous<Layout,Stride>, Layout> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
+  struct DynRankDimTraits<bool> {
 #endif
-
-  template <class Layout, unsigned Stride>
-  struct is_same< Layout, Kokkos::LayoutContiguous<Layout,Stride> > {
-    static const bool value = true;
+    using drdtraits = DynRankDimTraits<void>;
+    enum : size_t { unspecified = drdtraits::unspecified };
+
+    KOKKOS_INLINE_FUNCTION
+    static size_t computeRank(const size_t N0, const size_t N1, const size_t N2,
+                              const size_t N3, const size_t N4, const size_t N5,
+                              const size_t N6, const size_t N7) {
+      return drdtraits::computeRank(N0,N1,N2,N3,N4,N5,N6,N7);
+    }
+
+    template <typename Layout>
+    KOKKOS_INLINE_FUNCTION static size_t computeRank(const Layout& layout) {
+      return drdtraits::computeRank(layout);
+    }
+
+    template <typename Layout, typename... P>
+    KOKKOS_INLINE_FUNCTION static size_t computeRank(
+        const Kokkos::Impl::ViewCtorProp<P...>& prop ,
+        const Layout& layout) {
+      return drdtraits::computeRank(prop, layout);
+    }
+
+    // Non-contiguous Layout
+    template <typename Layout>
+    KOKKOS_INLINE_FUNCTION static std::enable_if_t<
+        !(is_layout_contiguous<Layout>::value),
+        Layout>
+    createLayout(const Layout& layout,
+                size_t new_rank = unspecified) {
+      return drdtraits::createLayout(layout, new_rank);
+    }
+
+    // Contiguous Layout
+    template <typename Layout>
+    KOKKOS_INLINE_FUNCTION static std::enable_if_t<
+        (is_layout_contiguous<Layout>::value),
+        Layout >
+    createLayout(const Layout& layout,
+                size_t new_rank = unspecified) {
+      return Layout(drdtraits::createLayout(layout.base_layout(), new_rank));
+    }
+
+    template <typename Traits, typename... P>
+    KOKKOS_INLINE_FUNCTION static std::enable_if_t<
+        !(is_layout_contiguous<typename Traits::array_layout>::value),
+        typename Traits::array_layout>
+    createLayout(const Kokkos::Impl::ViewCtorProp<P...>& prop,
+                typename Traits::array_layout layout) {
+#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
+      if constexpr (Traits::impl_is_customized &&
+                    !Kokkos::Impl::ViewCtorProp<P...>::has_accessor_arg) {
+        auto rank              = computeRank(prop, layout) - 1;
+        layout.dimension[rank] = unspecified;
+      }
+#endif
+      return createLayout(layout);
+    }
+
+    template <typename Traits, typename... P>
+    KOKKOS_INLINE_FUNCTION static std::enable_if_t<
+        (is_layout_contiguous<typename Traits::array_layout>::value),
+        typename Traits::array_layout>
+    createLayout(const Kokkos::Impl::ViewCtorProp<P...>& prop,
+                typename Traits::array_layout layout) {
+#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
+      if constexpr (Traits::impl_is_customized &&
+                    !Kokkos::Impl::ViewCtorProp<P...>::has_accessor_arg) {
+        auto rank              = computeRank(prop, layout) - 1;
+        layout.dimension[rank] = unspecified;
+      }
+#endif
+      return createLayout(layout);
+    }
+
+    template <typename ViewType, typename ViewArg>
+    static ViewType createView(const ViewArg& arg, const size_t N0,
+                              const size_t N1, const size_t N2, const size_t N3,
+                              const size_t N4, const size_t N5, const size_t N6,
+                              const size_t N7) {
+      return drdtraits::createView(arg, N0, N1, N2, N3, N4, N5, N6, N7);
+    }
   };
 
-  template <class Layout, unsigned Stride>
-#if defined(KOKKOS_COMPILER_INTEL)
-  inline constexpr bool is_same_v< Layout, Kokkos::LayoutContiguous<Layout,Stride>> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
-#else
-  static constexpr bool is_same_v< Layout, Kokkos::LayoutContiguous<Layout,Stride>> = is_same<Kokkos::LayoutContiguous<Layout,Stride>, Layout>::value;
-#endif
-}
+  // Overload reconstructLayout for LayoutContiguous
+  template <typename Layout, unsigned Stride, typename iType>
+  KOKKOS_INLINE_FUNCTION LayoutContiguous<Layout,Stride>
+  reconstructLayout(const LayoutContiguous<Layout,Stride>& layout, iType dynrank) {
+    return LayoutContiguous<Layout,Stride>(reconstructLayout(layout.base_layout(), dynrank));
+  }
+
+} // namespace Impl
+
+} // namespace Kokkos
 
 #if KOKKOS_VERSION >= 40499
 #include "View/Kokkos_ViewMapping.hpp"
diff --git a/packages/sacado/src/Kokkos_LayoutNatural.hpp b/packages/sacado/src/Kokkos_LayoutNatural.hpp
index 1a5ae982295..4cdd3228f99 100644
--- a/packages/sacado/src/Kokkos_LayoutNatural.hpp
+++ b/packages/sacado/src/Kokkos_LayoutNatural.hpp
@@ -64,21 +64,6 @@ struct inner_layout< LayoutNatural<Layout> > {
 
 } // namespace Kokkos
 
-// Make LayoutNatural<Layout> equivalent to Layout
-namespace std {
-
-  template <class Layout>
-  struct is_same< Kokkos::LayoutNatural<Layout>, Layout> {
-    static const bool value = true;
-  };
-
-  template <class Layout>
-  struct is_same< Layout, Kokkos::LayoutNatural<Layout> > {
-    static const bool value = true;
-  };
-
-}
-
 #if KOKKOS_VERSION >= 40499
 #include "View/Kokkos_ViewMapping.hpp"
 #else
diff --git a/packages/sacado/src/Sacado_Fad_Kokkos_Specialization.hpp b/packages/sacado/src/Sacado_Fad_Kokkos_Specialization.hpp
index eaf27171646..71877bf75ed 100644
--- a/packages/sacado/src/Sacado_Fad_Kokkos_Specialization.hpp
+++ b/packages/sacado/src/Sacado_Fad_Kokkos_Specialization.hpp
@@ -148,6 +148,33 @@ subview(const DynRankView<D, Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>,
   }
 }
 
+// Helper function to create a mapping for subdynrankview that avoids calling SubT::layout()
+// (which doesn't work correctly for LayoutContiguous<LayoutStride> and is slated for removal)
+template<class MapT, class SubT, size_t ... Idx>
+KOKKOS_INLINE_FUNCTION MapT
+create_subdynrankview_mapping(const SubT& sub, std::index_sequence<Idx...>) {
+  using return_layout_t = typename MapT::layout_type;
+  using exts_t = typename MapT::extents_type;
+  using idx_t = typename MapT::index_type;
+  const exts_t extents((Idx<SubT::rank() ? sub.extent(Idx):1)...);
+  if constexpr (std::is_same_v<return_layout_t, layout_stride>)
+  {
+    idx_t strides[7] = {(Idx<SubT::rank() ? sub.stride(Idx):1)...};
+    return MapT(mdspan_non_standard_tag(), extents, strides);
+  } else if constexpr(std::is_same_v<return_layout_t, Kokkos::Experimental::layout_left_padded<Kokkos::dynamic_extent> > ||
+                      std::is_same_v<return_layout_t, Kokkos::Experimental::layout_right_padded<Kokkos::dynamic_extent> >)
+  {
+    idx_t stride = 1;
+    if constexpr (SubT::rank() > 1) {
+      if constexpr(std::is_same_v<return_layout_t, Kokkos::Experimental::layout_left_padded<Kokkos::dynamic_extent> >)
+        stride = sub.stride(1);
+      else
+        stride = sub.stride(SubT::rank()-2);
+    }
+    return MapT(extents, stride);
+  }
+}
+
 template <class T, class LayoutSrc, unsigned StrideSrc, class... DRVArgs,
           class SubArg0 = int, class SubArg1 = int, class SubArg2 = int,
           class SubArg3 = int, class SubArg4 = int, class SubArg5 = int,
@@ -175,18 +202,11 @@ subdynrankview(const DynRankView<T, LayoutContiguous<LayoutSrc, StrideSrc>,
                                     // StrideSrc>,
       typename sub_t::device_type, typename sub_t::memory_traits>;
 
-  auto layout = sub.layout().base_layout();
-  for (int i = new_rank; i < 8; i++)
-    layout.dimension[i] = 1;
-  if constexpr (std::is_same_v<decltype(layout), LayoutStride>)
-    for (int i = new_rank; i < 8; i++)
-      layout.stride[i] = 1;
-
+  using return_map_t = typename return_type::mapping_type;
   return return_type{
       typename return_type::view_type(
           sub.data_handle(),
-          Impl::mapping_from_array_layout<typename return_type::mapping_type>(
-              layout),
+          create_subdynrankview_mapping<return_map_t>(sub, std::make_index_sequence<7>()),
           sub.accessor()),
       new_rank};
 }
diff --git a/packages/sacado/test/UnitTests/Fad_KokkosTests.hpp b/packages/sacado/test/UnitTests/Fad_KokkosTests.hpp
index e36c0f25310..da54757e38a 100644
--- a/packages/sacado/test/UnitTests/Fad_KokkosTests.hpp
+++ b/packages/sacado/test/UnitTests/Fad_KokkosTests.hpp
@@ -1430,7 +1430,8 @@ TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
 TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
   Kokkos_View_Fad, DynRankDimensionScalar, FadType, Layout, Device )
 {
-  typedef Kokkos::DynRankView<double,Layout,Device> DoubleViewType;
+  typedef typename Kokkos::inner_layout<Layout>::type DoubleLayout; // extract inner layout from LayoutContiguous
+  typedef Kokkos::DynRankView<double,DoubleLayout,Device> DoubleViewType;
   typedef Kokkos::DynRankView<FadType,Layout,Device> FadViewType;
   typedef typename FadViewType::size_type size_type;
"""


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "trilinos"
        self.version = "16.2.0"
        self.sha256 = "543aa56232d7c0cbe73705fab2d3b5524f11b15fef8917aa14de02d23a5ca418"
        self.filename = "trilinos-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-"
            + self.version.replace(".", "-")
            + ".tar.gz"
        )
        self.libraries = ["amesos2", "belos", "aztecoo", "amesos"]
        self.includes = [
            "Amesos2.hpp",
            "BelosSolverFactory.hpp",
            "Sacado.hpp",
            "AztecOO.h",
            "Amesos.h",
        ]
        self.dependencies = [
            "cmake",
            "openmpi",
            "lapack",
            "suitesparse",
            "superlu_dist",
            "mumps",
        ]

    def setDependencies(self, builder):
        return

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        if sys.platform == "darwin":
            with open(
                os.path.join(
                    builder._extract_dir, builder._extracted_folder, "sacado.patch"
                ),
                "w",
            ) as f:
                f.write(sacado_patch)

            builder.run_command(["patch", "-p1", "-i", "sacado.patch"])

        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
            builder.add_option("-DCMAKE_POSITION_INDEPENDENT_CODE=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")
        CC = builder.env["CC"]
        CXX = builder.env["CXX"]
        FC = builder.env["FC"]
        builder.add_option("-DCMAKE_C_COMPILER=" + CC)
        builder.add_option("-DCMAKE_CXX_COMPILER=" + CXX)
        builder.add_option("-DCMAKE_Fortran_COMPILER=" + FC)
        builder.add_option("-DTrilinos_SHOW_DEPRECATED_WARNINGS=OFF")
        builder.add_option("-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE")
        builder.add_option("-DTPL_ENABLE_Boost:BOOL=OFF")
        builder.add_option("-DTrilinos_SHOW_DEPRECATED_WARNINGS:BOOL=OFF")
        builder.add_option("-DTrilinos_ENABLE_SECONDARY_TESTED_CODE=ON")
        builder.add_option("-DTrilinos_ENABLE_Triutils:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_SEACAS:BOOL=OFF")
        builder.add_option("-DTrilinos_ENABLE_Epetra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Xpetra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Ifpack:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Ifpack2:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Teuchos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_ML:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_MueLu:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Stratimikos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Teko:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Belos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Amesos2:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Amesos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_AztecOO:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Sacado:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_EpetraExt:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Thyra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_ThyraEpetraAdapters:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Tpetra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Stratimikos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_TESTS:BOOL=OFF")
        builder.add_option("-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_MPI:BOOL=ON ")
        builder.add_option("-DMPI_BASE_DIR:PATH=" + builder.env["MPI_HOME"])
        builder.add_option("-DEpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_LAPACK:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_BLAS:BOOL=ON ")
        builder.add_option("-DHAVE_EPETRA_LAPACK_GSSVD3:BOOL=ON ")
        if "OPENBLAS_DIR" in builder.env:
            builder.add_option("-DLAPACK_LIBRARY_DIRS=" + builder.env["OPENBLAS_DIR"])
            builder.add_option("-DLAPACK_LIBRARY_NAMES=openblas")
            builder.add_option("-DBLAS_LIBRARY_DIRS=" + builder.env["OPENBLAS_DIR"])
            builder.add_option("-DBLAS_LIBRARY_NAMES=openblas")
        else:
            builder.add_option(
                "-DLAPACK_LIBRARY_DIRS=" + builder.env["LAPACK_DIR"] + "/lib"
            )
            builder.add_option("-DLAPACK_LIBRARY_NAMES=lapack;blas")
            builder.add_option(
                "-DBLAS_LIBRARY_DIRS=" + builder.env["LAPACK_DIR"] + "/lib"
            )
            builder.add_option("-DBLAS_LIBRARY_NAMES=blas")
        builder.add_option("-DTPL_ENABLE_UMFPACK:BOOL=ON ")
        builder.add_option(
            "-DUMFPACK_LIBRARY_NAMES:STRING=umfpack;amd;suitesparseconfig;cholmod;colamd;ccolamd;camd"
        )
        builder.add_option(
            "-DUMFPACK_LIBRARY_DIRS:PATH=" + builder.env["SUITESPARSE_DIR"] + "/lib"
        )
        builder.add_option(
            "-DUMFPACK_INCLUDE_DIRS:PATH="
            + builder.env["SUITESPARSE_DIR"]
            + "/include/suitesparse"
        )
        builder.add_option("-DTPL_ENABLE_AMD:BOOL=ON ")
        builder.add_option("-DAMD_LIBRARY_NAMES:STRING=amd;suitesparseconfig")
        builder.add_option(
            "-DAMD_LIBRARY_DIRS:PATH=" + builder.env["SUITESPARSE_DIR"] + "/lib"
        )
        builder.add_option(
            "-DAMD_INCLUDE_DIRS:PATH="
            + builder.env["SUITESPARSE_DIR"]
            + "/include/suitesparse"
        )
        if "SUPERLU_DIST_DIR" in builder.env:
            builder.add_option("-D Amesos_ENABLE_SuperLUDist:BOOL=ON ")
            builder.add_option("-DTPL_ENABLE_SuperLUDist:BOOL=ON ")
            builder.add_option("-DSuperLUDist_LIBRARY_NAMES:STRING=superlu_dist")
            builder.add_option(
                "-DSuperLUDist_LIBRARY_DIRS:PATH="
                + builder.env["SUPERLU_DIST_DIR"]
                + "/lib"
            )
            builder.add_option(
                "-DSuperLUDist_INCLUDE_DIRS:PATH="
                + builder.env["SUPERLU_DIST_DIR"]
                + "/include"
            )
        ext = utils.get_library_extension(builder.build_shared)
        if "PARMETIS_DIR" in builder.env:
            builder.add_option("-DTPL_ENABLE_ParMETIS:BOOL=ON ")
            builder.add_option("-D Amesos_ENABLE_ParMETIS:BOOL=ON ")
            builder.add_option(
                "-DParMETIS_LIBRARY_DIRS:PATH="
                + builder.env["PARMETIS_DIR"]
                + "/lib;"
                + builder.env["METIS_DIR"]
                + "/lib"
            )
            builder.add_option(
                "-DTPL_ParMETIS_INCLUDE_DIRS:PATH="
                + builder.env["PARMETIS_DIR"]
                + "/include;"
                + builder.env["METIS_DIR"]
                + "/include"
            )
            builder.add_option(
                "-DTPL_ParMETIS_LIBRARIES="
                + builder.env["PARMETIS_DIR"]
                + "/lib/libparmetis"
                + ext
                + ";"
                + builder.env["METIS_DIR"]
                + "/lib/libmetis"
                + ext
            )
        else:
            builder.add_option("-DTPL_ENABLE_ParMETIS:BOOL=OFF")
        builder.add_option("-DTPL_ENABLE_MUMPS:BOOL=ON ")
        builder.add_option(
            "-DMUMPS_LIBRARY_NAMES:STRING=dmumps;mumps_common;pord;scalapack;ptesmumps;ptscotch;ptscotcherr;scotch;scotcherr"
        )
        builder.add_option(
            "-DMUMPS_LIBRARY_DIRS:PATH="
            + builder.env["MUMPS_DIR"]
            + "/lib;"
            + builder.env["SCALAPACK_DIR"]
            + "/lib;"
            + builder.env["SCOTCH_DIR"]
            + "/lib"
        )
        builder.add_option(
            "-DMUMPS_INCLUDE_DIRS:PATH=" + builder.env["MUMPS_DIR"] + "/include"
        )
        builder.add_option("-DCMAKE_CXX_FLAGS:STRING=-DMUMPS_5_0")
        builder.add_option("-DAmesos_ENABLE_SCALAPACK:BOOL=ON ")
        builder.add_option(
            "-DSCALAPACK_LIBRARY_DIRS:FILEPATH=" + builder.env["SCALAPACK_DIR"] + "/lib"
        )
        builder.add_option("-DSCALAPACK_LIBRARY_NAMES:STRING=scalapack")
        builder.add_option("-D Amesos_ENABLE_LAPACK:BOOL=ON ")
        builder.add_option("-D Amesos_ENABLE_KLU:BOOL=ON ")
        builder.add_option("-D Amesos_ENABLE_UMFPACK:BOOL=ON ")
        builder.add_option("-D Amesos_ENABLE_MUMPS:BOOL=ON ")
        builder.add_option("-D Tpetra_INST_INT_INT:BOOL=ON ")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("TRILINOS_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
