from pyslimmc.operations import analysis_operation


def test_operation_descriptor_has_help():
    class X:
        @analysis_operation('demo')
        def f(self, x=1):
            return x
    x=X()
    assert x.f()==1
    assert x.f.help()=='demo'


def test_public_branch_classes_have_help_and_info():
    import pyslimmc._storage as s
    import pyslimmc.storage_analysis as a
    import pyslimmc.runs as r
    classes=[s.StorageRun,s.StorageSnapshots,s.StorageSnapshot,s.StorageStateSeries,
             s.StorageChains,s.StorageMomentsSeries,s.StorageChannelsSeries,
             s.StorageKineticsSeries,s.StorageActions,a.StorageFirings,
             a.StorageCopolymerization,a.StorageMicrostructure,s.StorageValidation,
             s.StorageDiagnostics,s.StorageRaw,r.Runs]
    for cls in classes:
        assert hasattr(cls,'help'), cls
        assert hasattr(cls,'info'), cls


def test_major_operations_have_pre_use_help():
    import pyslimmc._storage as s
    import pyslimmc.storage_analysis as a
    import pyslimmc.runs as r
    checks=[
      (s.StorageRun,['mwd','cld','chain_mass_spectrum','validate','mass_audit','summary']),
      (a.StorageCopolymerization,['mayo_lewis','reactivity_ratios','penultimate_parameters']),
      (a.StorageMicrostructure,['dyads','triads','blockiness']),
      (a.StorageFirings,['fire_shares','rate_shares']),
      (r.Runs,['filter','match','sweep','model_diff']),
    ]
    for cls,names in checks:
        for name in names:
            assert hasattr(getattr(cls,name),'help'), (cls,name)
