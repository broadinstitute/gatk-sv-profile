from __future__ import annotations

from gatk_sv_profile.aggregate import AggregatedData
from gatk_sv_profile.config import AnalysisConfig
from gatk_sv_profile.modules.base import AnalysisModule


class DemoModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "demo"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        del data, config


def test_analysis_module_output_dir_creates_directory(tmp_path) -> None:
    module = DemoModule()
    config = AnalysisConfig(output_dir=tmp_path)

    output_dir = module.output_dir(config)

    assert output_dir.exists()
    assert output_dir.name == "demo"
